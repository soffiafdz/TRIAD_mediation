#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lme4)
library(ADNIMERGE)
library(stargazer)
library(stringr)
#library(parameters)

ReFitModels <- F
USE_IMPUTED <- T

### FUNCTIONS
here("code/functions.R") |> source()

### INPUT
fpaths      <- here("data/rds", c("adni_hc-hvr.rds",
                                  "adni_dxs_imputed.rds",
                                  "adni_age_calculated.rds",
                                  "adni_cog-latent.rds",
                                  "adni_raket_edt.rds"))

## HCvol & HVR
if (file.exists(fpaths[1])) {
  hc_hvr.dt <- readRDS(fpaths[1])
} else {
  here("code/calc_hvr_adni.R") |> source()
}

## Imputed Dxs from ADNIMERGE
if (file.exists(fpaths[2])) {
  dx.dt     <- readRDS(fpaths[2])
} else {
  here("code/impute_dx.R") |> source()
}

## AGE
## Calculated age using the date of MRI sessions
if (file.exists(fpaths[3])) {
  age.dt    <- readRDS(fpaths[3])
} else {
  here("code/calculate_age_adni.R") |> source()
}

## Cognition
# Latent variable obtained from ADASQ4, MMSE, CDRSB and validated with CFA
if (file.exists(fpaths[4])) {
  cog.dt    <- readRDS(fpaths[4])
} else {
  here("code/cfa-cog_adni.R") |> source()
}

## RAKET EDT
if (file.exists(fpaths[5])) {
  raket.dt  <- readRDS(fpaths[5])
} else {
  here("code/parse_raket_adni.R") |> source()
}

rm(fpaths)


### Data PROCESSING
## Use imputed Dx?
if (USE_IMPUTED) dx.dt[, DX := stringr::str_remove(DX, "\\?")]
sel_dx.dt   <- dx.dt[, .(DX = factor(DX, levels = c("CN", "MCI", "Dementia"))),
                     .(PTID, EXAMDATE)]
#rm(USE_IMPUTED, dx.dt)

## PTGENDER & APOE4
data(adnimerge)
setDT(adnimerge)
tiv.dt      <- adnimerge[!duplicated(PTID), .(SEX = PTGENDER, APOE4), PTID]

## MERGE and average by side
DT          <- tiv.dt[age.dt, on = "PTID"
                      ][sel_dx.dt, on = .(PTID, EXAMDATE)
                      ][hc_hvr.dt, on = .(PTID, EXAMDATE)
                      ][, .(HCv = mean(HCvol_adj), HVR = mean(HVR)),
                      .(PTID, EXAMDATE, DX, AGE, SEX, APOE4)] |>
na.omit() |>
setkey(PTID)

DT          <- DTclean(DT, scalevars = c("HCv", "HVR"),
                       ordervars = c("APOE4", "DX"), centervars = "AGE",
                       reference_controls = TRUE)

# Get only AGE.c at baseline
DT          <- DT[.(1), on = "VIS", .(PTID, AGE.bl.c = AGE.c)
                  ][DT, on = "PTID"
                  ] |> setcolorder("AGE.bl.c", after = "AGE.bl")

# Melting
DT_l        <- DT |>
melt(measure = patterns("^H(Cv|VR).scl$"),
     variable.name = "HC", value.name = "VAL")
hcvars      <- DT_l[, levels(HC)]

# Rename vars/covars for easier formatting:
vars        <- c("VAL", "TIME", "SEX", "DX", "APOE4", "AGE.bl.c",
                 "VIS", "PTID")
vars_short  <- c("Y", "T", "S", "D", "A4", "A", "I", "ID")
setnames(DT_l, vars, vars_short)

### MLM with lmer: HC & HVR
## Y: HC integrity with HCvol & HVR
# Params and model definitions
# After iteration : skip 4, 9, 12, 14, 15
f1.1        <- build_formulas("Y", c("T", "D", "S", "A4", "A"),
                              random_effects = "slopes", slopevars = "I",
                              quadraticvars = "T", quadratic_interactions = T)

f1.2        <- build_formulas("Y", c("T", "D", "S", "A4", "A"),
                              random_effects = "slopes", slopevars = "I",
                              quadraticvars = "T", quadratic_interactions = T,
                              skip_items = c(8:11, 13:14))

# Fits
f1 <- f1.2
fpath <- here("data/rds/adni_mlm-hc-hvr.rds")
if (!file.exists(fpath) | ReFitModels) {
  fits1        <- vector("list", 2)
  names(fits1) <- hcvars
  for (i in 1:2) {
    fits1[[i]] <- vector("list", length(f1))
    for (j in seq_along(fits1[[i]])) {
      fits1[[i]][[j]] <- lmer(as.formula(f1[j]),
                             DT_l[hcvars[i], on = "HC"], REML = F,
                             control = lmerControl(optimizer = "bobyqa"))
    }
  }
  saveRDS(fits1, fpath)
} else {
  fits1        <- readRDS(fpath)
}

## LRT (anova)
# HC volume
# Best model (10):
# Y ~ T+D+S+A4+A+T2+T:D+T:S+T:A4+T:A+D:S+D:A4+T2:D+T2:A+rs
# do.call(anova, fits1$HCv.scl) |> suppressWarnings()

# HVR
# Best model (10):
# Y ~ T+D+S+A4+A+T2+T:D+T:S+T:A4+T:A+D:S+D:A4+T2:D+T2:A+rs
# do.call(anova, fits1$HVR.scl) |> suppressWarnings()

pred_labels <- c("HC volume", "HVR")
covs_labels <- fits1[[1]][[10]] |> summary() |> coef() |> rownames() |>
str_remove_all("\\(|\\)|\\.L") |> str_replace_all("\\^2", ".Q") |>
str_replace_all("(.+):IT.Q", "IT.Q:\\1") |> str_replace_all("I?T", "Time") |>
str_replace_all("(:?)\\.Q(:?)", "\\1(Q)\\2") |>
str_replace_all("^([^:]*)\\(Q\\)(?!:)", "\\1 (Quad.)") |>
str_replace_all("SMale", "Sex(M)") |> str_replace_all("A(?!4)", "Age(bl)") |>
str_replace_all("D", "Dx") |> str_replace_all("^A4(\\s?)", "APOE4\\1") |>
str_replace_all("^([^:]*)\\(bl\\)(?!:)", "\\1 (baseline)") |>
str_replace_all("(^[^:]*)\\(M\\)", "\\1 (Male)")

stargazer(fits1[[1]][10], fits1[[2]][10],
          out = here("data/derivatives/adni_mlm-hc-hvr.html"),
          intercept.bottom = F, type = "html", single.row = T,
          title = "MLM (random slopes: VISIT|ID)", column.separate = c(1, 2),
          dep.var.labels.include = FALSE, column.labels = pred_labels,
          dep.var.caption = "Hippocampal integrity",
          covariate.labels = covs_labels)

### EDT (Raket's method)
## Data PROCESSING
setkey(raket.dt, PTID, EXAMDATE)
DT_r        <- DT[, .(SEX, APOE4, HCv, HVR), keyby = .(PTID, EXAMDATE)
                  ][raket.dt]

# Time invariant covariates
tics.cols   <- c("PTID", "SEX", "APOE4")
tic.dt      <- DT_r[, ..tics.cols] |>
na.omit() |>
unique() |>
setkey(PTID)
DT_r[, (tics.cols[-1]) := NULL]

# Baseline HC integrity
bl.cols     <- c("EXAMDATE", "HCv", "HVR")
hchvr.bl.dt <- DT_r[!is.na(HCv), .SD[which.min(AGE)],
                    "PTID", .SDcols = bl.cols ] |>
setnames(bl.cols, paste0(bl.cols, ".bl")) |>
setkey(PTID)

DT_r        <- tic.dt[hchvr.bl.dt
                      ][DT_r
                      ][EXAMDATE.bl <= EXAMDATE] |>
unique() |>
setcolorder(c("EXAMDATE", "MONTH", "AGE"), after = "HVR.bl")
DT_r[, EXAMDATE.bl := NULL]

# Final CLEAN
DT_r        <- DTclean(DT_r, scalevars = c("HCv", "HVR", "HCv.bl", "HVR.bl"),
                       centervars = "AGE")
DT_r        <- DT_r[VIS == 1, .(PTID, MONTH.bl = MONTH)][DT_r, on = "PTID"]
DT_r[, MONTH := MONTH - MONTH.bl]
DT_r[, MONTH.bl := NULL]

# Get only AGE.c at baseline
DT_r        <- DT_r[.(1), on = "VIS", .(PTID, AGE.bl.c = AGE.c)
                    ][DT_r, on = "PTID"
                    ] |> setcolorder("AGE.bl.c", after = "AGE.bl")

# Melting
DT_lr       <- DT_r |>
melt(measure = patterns(VAL = "(HCv|HVR).scl", VAL.bl = "(HCv|HVR).bl.scl"),
     variable.name = "HC")
DT_lr[, HC := factor(HC, labels = hcvars)]

# Rename vars/covars for easier formatting:
vars        <- c("EDT", "VAL", "VAL.bl", "AGE.bl.c", "SEX", "APOE4", "VIS", "PTID")
vars_short  <- c("Y", "X", "Xb", "A", "S", "A4", "I", "ID")
setnames(DT_lr, vars, vars_short)

## Params and model definitions
# 1:7 — Y: EDT (baseline HCintegrity)
f2          <- build_formulas("Y", c("Xb", "A", "S", "A4"),
                              random_effects = "slopes", slopevars = "I")
f2          <- c(f2, gsub("Xb", "X", f2)) # 8:14 — Y: EDT (long HCintegrity)

## Fits
fpath <- here("data/rds/adni_mlm-edt.rds")
if (!file.exists(fpath) | ReFitModels) {
  fits2        <- vector("list", 2)
  names(fits2) <- hcvars
  for (i in 1:2) {
    fits2[[i]] <- vector("list", length(f2))
    for (j in seq_along(f2)) {
      fits2[[i]][[j]] <- lmer(as.formula(f2[j]),
                             DT_lr[hcvars[i], on = "HC"], REML = F,
                             control = lmerControl(optimizer = "bobyqa"))
    }
  }
  saveRDS(fits2, fpath)
} else {
  fits2        <- readRDS(fpath)
}

## LRT (anova)
# HC volume
# Best model (4) — EDT (baseline HCintegrity)
# Y ~ X+A+S+A4+X:A+X:S+X:A4+ri
# Best model (11) — EDT (longitudinal HCintegrity)
# Y ~ X+A+S+A4+ri
 #do.call(anova, fits2$HCv.scl[1:7]) |> suppressWarnings()
 #do.call(anova, fits2$HCv.scl[8:14]) |> suppressWarnings()

# HVR
# Best model (1) — EDT (baseline HCintegrity)
# Y ~ X+A+S+A4+X:A+X:S+X:A4+ri
# Best model (9) — EDT (longitudinal HCintegrity)
# Y ~ X+A+S+A4+X:Ari
 #do.call(anova, fits2$HVR.scl[1:7]) |> suppressWarnings()
 #do.call(anova, fits2$HVR.scl[8:14]) |> suppressWarnings()

### Out table
pred_labels <- c("Est. Disease Time", "HC volume", "HVR")
covs_labels <- fits2[[1]][[4]] |> summary() |> coef() |> rownames() |>
str_remove_all("\\(|\\)|\\.L") |>
str_replace_all("(:?)\\.Q(:?)", "\\1(Q)\\2") |>
str_replace_all("^([^:]*)\\(Q\\)(?!:)", "\\1 (Quad.)") |>
str_replace_all("SMale", "Sex(M)") |> str_replace_all("A(?!4)", "Age(bl)") |>
str_replace_all("^A4(\\s?)", "APOE4\\1") |> str_replace_all("Xb", "HC(bl)") |>
str_replace_all("^([^:]*)\\(bl\\)(?!:)", "\\1 (baseline)") |>
str_replace_all("(^[^:]*)\\(M\\)", "\\1 (Male)")

stargazer(fits2[[1]][4], fits2[[2]][1], fits2[[2]][4],
          out = here("data/derivatives/adni_mlm-edt1.html"),
          intercept.bottom = F, type = "html", single.row = T,
          title = "MLM: (random slopes: VISIT|ID)", column.separate = c(1, 2),
          dep.var.labels.include = FALSE, column.labels = pred_labels[-1],
          dep.var.caption = pred_labels[1], covariate.labels = covs_labels)

covs_labels <- str_remove_all(covs_labels, "(?<=HC)\\s?\\(\\w*\\)")
stargazer(fits2[[1]][11], fits2[[2]][9], fits2[[2]][11],
          out = here("data/derivatives/adni_mlm-edt2.html"),
          intercept.bottom = F, type = "html", single.row = T,
          title = "MLM: (random slopes: VISIT|ID)", column.separate = c(1, 2),
          dep.var.labels.include = FALSE, column.labels = pred_labels[-1],
          dep.var.caption = pred_labels[1], covariate.labels = covs_labels)

### Predict COGNITION
cog.dt      <- age.dt[cog.dt[, .(PTID, VISCODE, COG_latent)],
                      on = .(PTID, VISCODE)
                      ][, COG_latent, keyby = .(PTID, EXAMDATE)]
DT_c        <- DT[, .(SEX, AGE, APOE4, HCv, HVR), keyby = .(PTID, EXAMDATE)
                  ][cog.dt] |> setkey(PTID)

# Time invariant covariates
tics.cols   <- c("PTID", "SEX", "APOE4")
tic.dt      <- DT_c[, ..tics.cols] |>
na.omit() |>
unique() |>
setkey(PTID)
DT_c[, (tics.cols[-1]) := NULL]

# Baseline HC integrity
hchvr.bl.dt <- DT_c[!is.na(HCv), .SD[which.min(AGE)],
                    "PTID", .SDcols = bl.cols ] |>
setnames(bl.cols, paste0(bl.cols, ".bl")) |>
setkey(PTID)

DT_c        <- tic.dt[hchvr.bl.dt][DT_c][EXAMDATE.bl <= EXAMDATE] |>
unique()
DT_c[, EXAMDATE.bl := NULL]

# Shift Cog factor to being positive
DT_c[, COG := COG_latent + abs(min(COG_latent))]

DT_c        <- DTclean(DT_c, ordervars = "APOE4",
                       scalevars = c("HCv", "HCv.bl", "HVR", "HVR.bl", "COG"),
                       centervars = "AGE")
# AGE.c at baseline
DT_c        <- DT_c[.(1), on = "VIS", .(PTID, AGE.bl.c = AGE.c)
                    ][DT_c, on = "PTID"
                    ] |> setcolorder("AGE.bl.c", after = "AGE.bl")

# Melting
DT_lc       <- DT_c |>
melt(measure = patterns(VAL = "(HCv|HVR).scl", VAL.bl = ".bl.scl"),
     variable.name = "HC")
DT_lc[, HC := factor(HC, labels = hcvars)]

# Rename vars/covars for easier formatting:
vars        <- c("COG.scl", "VAL", "VAL.bl",
                 "TIME", "SEX", "APOE4", "AGE.bl.c",
                 "VIS", "PTID")
vars_short  <- c("Y", "X", "Xb", "T", "S", "A4", "A", "I", "ID")
setnames(DT_lc, vars, vars_short)

## Models
# 1:16  — Y: COG (baseline HCintegrity)
f3.1        <- build_formulas("Y", c("T", "Xb", "S", "A4", "A"),
                            random_effects = "slopes", slopevars = "I",
                            quadraticvars = "T", quadratic_interactions = T)
f3.2        <- build_formulas("Y", c("T", "Xb", "S", "A4", "A"),
                            random_effects = "slopes", slopevars = "I",
                            quadraticvars = "T", quadratic_interactions = T,
                            skip_items = c(9:11, 13))

# 17:32 — Y: COG (baseline & longitudinal HCintegrity)
f3.1        <- c(f3.1, gsub("Xb", "X", f3.1))
f3.2        <- c(f3.2, gsub("Xb", "X", f3.2))

f3 <- f3.2

## Fits
fpath       <- here("data/rds/adni_mlm-cog.rds")
if (!file.exists(fpath) | ReFitModels) {
  fits3         <- vector("list", 2)
  names(fits3)  <- hcvars ## Defined above
  for (i in 1:2) {
    fits3[[i]]  <- vector("list", length(f3))
    for (j in seq_along(f3)) {
      fits3[[i]][[j]] <- lmer(as.formula(f3[[j]]),
                              DT_lc[hcvars[i], on = "HC"], REML = F,
                             control = lmerControl(optimizer = "bobyqa"))
    }
  }
  saveRDS(fits3, fpath)
} else {
  fits3         <- readRDS(fpath)
}

## LRT (anova)
# HC volume
# Best model (12) — COG (baseline HCintegrity)
# Y ~ T+X+S+A4+A+T2+T:X+T:S+T+A4+T:A+X:S+X:A4+T2:X+T2:A4+T2:A+rs
# Best model (22) — COG (baseline & longitudinal HCintegrity)
# Y ~ T+X+S+A4+A+T2+T:X+T:S+T+A4+T:A+X:S+X:A4+T2:X+T2:A4+T2:A+rs
#do.call(anova, fits3$HCv.scl[1:16]) |> suppressWarnings()
#do.call(anova, fits3$HCv.scl[17:32]) |> suppressWarnings()

# HVR
# Best model (12) — COG (baseline HCintegrity)
# Y ~ T+X+S+A4+A+T2+T:X+T:S+T+A4+T:A+X:S+X:A4+T2:X+rs
# Best model (24) — COG (baseline & longitudinal HCintegrity)
# Y ~ T+X+S+A4+A+T2+T:X+T:S+T+A4+T:A+X:S+X:A4+T2:X+T2+A4+T2:A+rs
#do.call(anova, fits3$HVR.scl[1:16]) |> suppressWarnings()
#do.call(anova, fits3$HVR.scl[17:32]) |> suppressWarnings()

pred_labels <- c("Cognition (latent variable)", "HC volume", "HVR")
covs_labels <- fits3[[1]][[12]] |> summary() |> coef() |> rownames() |>
str_remove_all("\\(|\\)|\\.L") |> str_replace_all("\\^2", ".Q") |>
str_replace_all("(.+):IT.Q", "IT.Q:\\1") |> str_replace_all("I?T", "Time") |>
str_replace_all("(:?)\\.Q(:?)", "\\1(Q)\\2") |>
str_replace_all("^([^:]*)\\(Q\\)(?!:)", "\\1 (Quad.)") |>
str_replace_all("Xb", "HC(bl)") |> str_replace_all("SMale", "Sex(M)") |>
str_replace_all("A(?!4)", "Age(bl)") |>
str_replace_all("D", "Dx") |> str_replace_all("^A4(\\s?)", "APOE4\\1") |>
str_replace_all("^([^:]*)\\(bl\\)(?!:)", "\\1 (baseline)") |>
str_replace_all("(^[^:]*)\\(M\\)", "\\1 (Male)")

stargazer(fits3[[1]][12], fits3[[2]][12],
          out = here("data/derivatives/adni_mlm-cog1.html"),
          intercept.bottom = F, type = "html", single.row = T,
          title = "MLM (random slopes: VISIT|ID)",
          dep.var.labels.include = FALSE, column.labels = pred_labels[-1],
          dep.var.caption = pred_labels[1], covariate.labels = covs_labels)

covs_labels <- str_remove_all(covs_labels, "(?<=HC)\\s?\\(\\w*\\)") |>
str_replace_all("HC", "HVR")
#stargazer(fits3[[1]][22], fits3[[1]][24], fits3[[2]][24],
stargazer(fits3[[2]][24],
          out = here("data/derivatives/adni_mlm-cog2.html"),
          intercept.bottom = F, type = "html", single.row = T,
          title = "MLM (random slopes: VISIT|ID)",
          dep.var.labels.include = FALSE,
          #column.labels = pred_labels[-1],
          #column.separate = c(2, 1),
          dep.var.caption = pred_labels[1], covariate.labels = covs_labels)
