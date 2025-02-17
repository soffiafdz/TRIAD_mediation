#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lme4)
library(ADNIMERGE)
library(stringr)
library(stargazer)
library(progress)
#library(parameters)

REFITMODELS <- F
PRINTPLOTS <- T

### FUNCTIONS
here('code/functions.R') |> source()

### INPUT
fpaths <- list(
  RDS = c(
  "adni_hc-hvr",
  "adni_dxs_imputed",
  "adni_age_calculated",
  "ucb_pet-amy",
  "ucb_pet-tau",
  "../adni_wmh-vol",
  "adni_cog-latent",
  "adni_biomarkers-latent"
  ) |> sprintf(fmt = "data/rds/%s.rds") |> here(),
  SRC = c(
    "calc_hvr_adni",
    "impute_dx",
    "calculate_age_adni",
    "parse_pet_adni",
    "cfa-cog_adni",
    "cfa-biomarkers_adni"
  ) |> sprintf(fmt = "code/%s.R") |> here()
)

## HCvol & HVR
if (file.exists(fpaths[["RDS"]][1])) {
  hc_hvr.dt <- fpaths[["RDS"]][1] |> readRDS()
} else {
  fpaths[["SRC"]][1] |> source()
}

## Imputed Dxs from ADNIMERGE
if (file.exists(fpaths[["RDS"]][2])) {
  dx.dt     <- fpaths[["RDS"]][2] |> readRDS()
} else {
  fpaths[["SRC"]][2] |> source()
}

## AGE
## Calculated age using the date of MRI sessions
if (file.exists(fpaths[["RDS"]][3])) {
  age.dt    <- fpaths[["RDS"]][3] |> readRDS()
} else {
  fpaths[["SRC"]][3] |> source()
}

## PET
if (all(file.exists(fpaths[["RDS"]][4:5]))) {
  ucb_amy.dt <- fpaths[["RDS"]][4] |> readRDS()
  ucb_tau.dt <- fpaths[["RDS"]][5] |> readRDS()
} else {
  fpaths[["SRC"]][4] |> source()
}

## WMH
if (
  fpaths[["RDS"]][6] |>
    sub(pattern = "\\.rds", replacement = "\\.csv") |>
    file.exists()
) {
  # TODO: Check that gsub works
  wmh.dt <- fpaths[["RDS"]][6] |>
    sub(pattern = "\\.rds", replacement = "\\.csv") |>
    fread()
} else {
  sprintf("File: %s is required but could not be found.", fpaths[6]) |> stop()
}

## Cognition
# Latent variable obtained from ADASQ4, MMSE, CDRSB and validated with CFA
if (file.exists(fpaths[["RDS"]][7])) {
  cog.dt    <- fpaths[["RDS"]][7] |> readRDS()
} else {
  fpaths[["SRC"]][5] |> source()
}

## Pathology
# Latent variables obtained from UCBerkeley Amy, Tau PET and WMH vols
if (file.exists(fpaths[["RDS"]][8])) {
  path.dt   <- fpaths[["RDS"]][8] |> readRDS()
} else {
  fpaths[["SRC"]][6] |> source()
}
rm(fpaths)

### Data PROCESSING
## PTGENDER & APOE4
data(adnimerge)
setDT(adnimerge)
tiv.dt <- adnimerge[!duplicated(PTID), .(SEX = PTGENDER, APOE4), PTID]


## Base Data.table
# MERGE and average by side
DT <- tiv.dt[
  age.dt,
  on = "PTID"
][
  hc_hvr.dt,
  on = .(PTID, EXAMDATE)
][
  ,
  .(HCv = mean(HCvol_adj), HVR = mean(HVR)),
  .(PTID, VISCODE, EXAMDATE, AGE, SEX, APOE4)
][
  wmh.dt,
  on = .(PTID, EXAMDATE)
  -c("VISCODE")
] |>
  na.omit() |>
  setkey(PTID)

### PET
## Amyloid of only subjects with MRI
ab.dt <- ucb_amy.dt[PTID %in% DT[!duplicated(PTID), PTID], c(1,3,7)] |>
  setnames("CENTILOIDS", "AMY")
ab.dt[, EXAMDATE := DATE_AMY]

## DT with only earliest session
#ab.early.dt <- ab.dt[
  #,
  #.SD[which.min(DATE_AMY)],
  #PTID,
  #.SDcols = c("DATE_AMY", "AMY")
#]

## Tau of subjects with MRI
tau.dt <- ucb_tau.dt[PTID %in% DT[!duplicated(PTID), PTID], c(1,3,7)] |>
  setnames("META_TEMPORAL_SUVR", "TAU")
tau.dt[, EXAMDATE := DATE_TAU]

## DTs with earliest and latest
#tau.early.dt<- tau.dt[, .SD[which.min(DATE_TAU)], PTID,
                      #.SDcols = c("DATE_TAU", "TAU")]
#tau.late.dt <- tau.dt[, .SD[which.max(DATE_TAU)], PTID,
                      #.SDcols = c("DATE_TAU", "TAU")]

### TODO: MOVE THIS TO ITS RESPECTIVE SECTION
### CSF biomarkers
#csf.cols    <- c("AMY_csf", "TAU_csf", "PTAU_csf")
#csf.dt      <- adnimerge[, .(ABETA, TAU, PTAU), .(PTID, EXAMDATE)
                         #][PTID %in% DT[!duplicated(PTID), PTID]] |>
#setnames(3:5, csf.cols) |> na.omit()

#csf.dt[AMY_csf %like% ">|<", AMY_csf := stringr::str_remove(AMY_csf, ">|<")]
#csf.dt[TAU_csf %like% ">|<", TAU_csf := stringr::str_remove(TAU_csf, ">|<")]
#csf.dt[PTAU_csf %like% ">|<", PTAU_csf := stringr::str_remove(PTAU_csf, ">|<")]
#csf.dt[, (csf.cols) := lapply(.SD, as.numeric), .SDcols = csf.cols]

DT_p <- ab.dt[DT, on = .(PTID, EXAMDATE), roll = T]
DT_p <- tau.dt[DT_p, on = .(PTID, EXAMDATE), roll = T] |>
  na.omit() |>
  setcolorder(c("DATE_AMY", "DATE_TAU"), after = "EXAMDATE") |>
  setcolorder(c("AMY", "TAU"), after = "WMHvol")

#rm(tiv.dt, hc_hvr.dt, wmh.dt, ucb_amy.dt, ucb_tau.dt)

## TODO: Move this to its part
##[path.dt, on = .(PTID, EXAMDATE)] |>

## Data cleaning
DT_p[, DATE_AMY := EXAMDATE - DATE_AMY]
DT_p[, DATE_TAU := EXAMDATE - DATE_TAU]
setnames(DT_p, c("DATE_AMY", "DATE_TAU"), c("AMY.diff", "TAU.diff"))
#DT_p[, DX := factor(DX, levels = c("CN", "MCI", "Dementia"))]
DT_p <- DTclean(
  DT_p,
  scalevars = names(DT_p)[8:12],
  ordervars = "APOE4",
  centervars = "AGE"
)

# Get only AGE.c at baseline
DT_p <- DT_p[
  .(1),
  on = "VIS",
  .(PTID, AGE.bl.c = AGE.c)
][
  DT_p,
  on = "PTID"
] |> setcolorder("AGE.bl.c", after = "AGE.bl")

# Melting
DT_lp <- melt(
  DT_p,
  measure = patterns("^H(Cv|VR).scl$"),
  variable.name = "HC",
  value.name = "VAL"
)
hcvars      <- DT_lp[, levels(HC)]

## Rename vars/covars for easier formatting:
vars <- c(
  "VAL", "TIME", "SEX", "APOE4", "AGE.bl.c",
  "WMHvol.scl", "AMY.scl", "AMY.diff", "VIS", "PTID"
)

vars_short <- c(
  "Y", "T", "S", "A4", "A", "W", "P", "Pt", "I", "ID"
)

setnames(DT_lp, vars, vars_short)


### MLM with lmer
## Extend HCv / HVR models with pathology
## Params and model definitions
f1 <- build_formulas(
  "Y",
  c("T", "S", "A4", "W", "A", "P", "Pt"),
  notinteractionvars = c("A", "Pt"),
  quadraticvars = "T",
  quadratic_interactions = T,
  random_effects = "intercepts",
  skip_items = c(4,6,11,14,15,17)
)

### Fits
fpath <- here("data/rds/adni_mlm-pet-amyloid.rds")
if (!file.exists(fpath) | ReFitModels) {
#if (TRUE) {
  fits1 <- vector("list", 2)
  names(fits1) <- hcvars
  for (i in 1:2) {
    fits1[[i]] <- vector("list", length(f1))
    for (j in seq_along(f1)) {
      fits1[[i]][[j]] <- lmer(
        as.formula(f1[j]),
        DT_lp[hcvars[i], on = "HC"],
        REML = F,
        control = lmerControl(optimizer = "bobyqa")
      )
    }
  }
  saveRDS(fits1, fpath)
} else {
  fits1 <- readRDS(fpath)
}

#### LRT
## Model 9 is significantly better
## HC volume
# Y ~ X+T:A4+T:P+S:W+A4:W+A4:P+T2:S+T2:A4+T2:P
##do.call(anova, fits1$HCv_scld) |> suppressWarnings()

## HVR
# Y ~ X+T:A4+T:P+S:W+A4:P+T2:A4+T2:P
##do.call(anova, fits1$HVR_scld) |> suppressWarnings()

pred_labels <- c("HC volume", "HVR")
covs_labels <- fits1[[1]][[12]] |> summary() |> coef() |> rownames() |>
str_remove_all("\\(|\\)|\\.L") |> str_replace_all("\\^2", ".Q") |>
str_replace_all("(.+):IT.Q", "IT.Q:\\1") |> str_replace_all("I?T", "Time") |>
str_replace_all("(:?)\\.Q(:?)", "\\1(Q)\\2") |>
str_replace_all("^([^:]*)\\(Q\\)(?!:)", "\\1 (Quad.)") |>
str_replace_all("SMale", "Sex(M)") |> str_replace_all("A(?!4)", "Age(bl)") |>
str_replace_all("^P$", "PET (Amyloid)") |>str_replace_all("(:?)P$", "\\1PET") |>
str_replace_all("^Pt$", "PET/MRI time offset") |>
str_replace_all("W", "WMH") |> str_replace_all("^A4(\\s?)", "APOE4\\1") |>
str_replace_all("^([^:]*)\\(bl\\)(?!:)", "\\1 (baseline)") |>
str_replace_all("(^[^:]*)\\(M\\)", "\\1 (Male)")

#stargazer(fits1[[1]][10], fits1[[2]][10],
          #out = here("data/derivatives/adni_mlm-hc-pet-amy.html"),
          #intercept.bottom = F, type = "html", single.row = T,
          #title = "MLM (random intercepts)",
          #dep.var.labels.include = FALSE, column.labels = pred_labels,
          #dep.var.caption = "Hippocampal integrity",
          #covariate.labels = covs_labels)


f2 <- build_formulas(
  "Y",
  c("T", "S", "A4", "W", "A", "P", "Pt"),
  notinteractionvars = c("A", "Pt"),
  quadraticvars = "T",
  quadratic_interactions = T,
  random_effects = "intercepts",
  skip_items = c(4,6,9,10,12:15,17)
)

## Tau PET
setnames(
  DT_lp,
  c("P", "Pt", "TAU.scl", "TAU.diff"),
  c("AMY.scl", "AMY.diff", "P", "Pt")
)

fpath <- here("data/rds/adni_mlm-pet-tau.rds")
if (!file.exists(fpath) | ReFitModels) {
  fits2 <- vector("list", 2)
  names(fits2) <- hcvars
  for (i in 1:2) {
    fits2[[i]] <- vector("list", length(f2))
    for (j in seq_along(f2)) {
      fits2[[i]][[j]] <- lmer(
        as.formula(f2[j]),
        DT_lp[hcvars[i], on = "HC"],
        REML = F,
        control = lmerControl(optimizer = "bobyqa")
      )
    }
  }
  saveRDS(fits2, fpath)
} else {
  fits2 <- readRDS(fpath)
}

### LRT
## HC volume
## Y ~ X+T:A4+T:P+S:W+W:P
#do.call(anova, fits2$HCv_scld) |> suppressWarnings()

## HVR
## Y ~ X+T:S+T:A4+T:P+S:W+S:P+W:P+T2:P
#do.call(anova, fits2$HVR_scld) |> suppressWarnings()

pred_labels <- c("HC volume", "HVR")
covs_labels <- fits2[[1]][[9]] |> summary() |> coef() |> rownames() |>
str_remove_all("\\(|\\)|\\.L") |> str_replace_all("\\^2", ".Q") |>
str_replace_all("(.+):IT.Q", "IT.Q:\\1") |> str_replace_all("I?T", "Time") |>
str_replace_all("(:?)\\.Q(:?)", "\\1(Q)\\2") |>
str_replace_all("^([^:]*)\\(Q\\)(?!:)", "\\1 (Quad.)") |>
str_replace_all("SMale", "Sex(M)") |> str_replace_all("A(?!4)", "Age(bl)") |>
str_replace_all("^P$", "PET (Tau)") |>str_replace_all("(:?)P$", "\\1PET") |>
str_replace_all("^Pt$", "PET/MRI time offset") |>
str_replace_all("W", "WMH") |> str_replace_all("^A4(\\s?)", "APOE4\\1") |>
str_replace_all("^([^:]*)\\(bl\\)(?!:)", "\\1 (baseline)") |>
str_replace_all("(^[^:]*)\\(M\\)", "\\1 (Male)")

#stargazer(fits2[[1]][9], fits2[[2]][8], fits2[[2]][9],
          #out = here("data/derivatives/adni_mlm-hc-pet-tau.html"),
          #intercept.bottom = F, type = "html", single.row = T,
          #title = "MLM (random intercepts)", column.separate = c(1, 2),
          #dep.var.labels.include = FALSE, column.labels = pred_labels,
          #dep.var.caption = "Hippocampal integrity",
          #covariate.labels = covs_labels)

# Both Tau & AMY
setnames(
  DT_lp,
  c("P", "Pt", "AMY.scl", "AMY.diff"),
  c("Ptau", "Pt1", "Pamy", "Pt2")
)

## Formula is manual
#Y ~
f3 <- paste0(
  "Y ~ T+S+A4+W+A+Ptau+Pamy+Pt1+Pt2+I(T^2)",
  "+T:S+T:A4+T:Pamy+T:Ptau+S:W+S:Pamy+S:Ptau+A4:Pamy+W:Ptau",
  "+I(T^2):S+I(T^2):A4+I(T^2):Pamy+I(T^2):Ptau+(1|ID)"
)

fpath <- here("data/rds/adni_mlm-pet.rds")
if (!file.exists(fpath) | ReFitModels) {
  fits3 <- vector("list", 2)
  names(fits3) <- hcvars
  for (i in 1:2) {
    fits3[[i]] <- vector("list", length(f3))
    fits3[[i]] <- lmer(
      as.formula(f3),
      DT_lp[hcvars[i], on = "HC"],
      REML = F,
      control = lmerControl(optimizer = "bobyqa")
    )
  }
  saveRDS(fits3, fpath)
} else {
  fits3 <- readRDS(fpath)
}

pred_labels <- c("HC volume", "HVR")
covs_labels <- fits3[[1]] |> summary() |> coef() |> rownames() |>
str_remove_all("\\(|\\)|\\.L") |> str_replace_all("\\^2", ".Q") |>
str_replace_all("(.+):IT.Q", "IT.Q:\\1") |> str_replace_all("I?T", "Time") |>
str_replace_all("(:?)\\.Q(:?)", "\\1(Q)\\2") |>
str_replace_all("^([^:]*)\\(Q\\)(?!:)", "\\1 (Quad.)") |>
str_replace_all("SMale", "Sex(M)") |> str_replace_all("A(?!4)", "Age(bl)") |>
str_replace_all("^Ptau$", "PET (Tau)") |>
str_replace_all("^Pamy$", "PET (Amy)") |>
str_replace_all("(:?)Ptau$", "\\1Tau") |>
str_replace_all("(:?)Pamy$", "\\1Amy") |>
str_replace_all("^Pt1$", "PET(Tau)/MRI offset") |>
str_replace_all("^Pt2$", "PET(Amy)/MRI offset") |>
str_replace_all("W", "WMH") |> str_replace_all("^A4(\\s?)", "APOE4\\1") |>
str_replace_all("^([^:]*)\\(bl\\)(?!:)", "\\1 (baseline)") |>
str_replace_all("^([^:]*)\\(M\\)(?!:)", "\\1 (Male)")

#stargazer(fits3[[1]], fits3[[2]],
          #out = here("data/derivatives/adni_mlm-hc-pet.html"),
          #intercept.bottom = F, type = "html", single.row = T,
          #title = "MLM (random intercepts)",
          #dep.var.labels.include = FALSE, column.labels = pred_labels,
          #dep.var.caption = "Hippocampal integrity",
          #covariate.labels = covs_labels)


## Predict COGNITION
cog.dt <- cog.dt[, .(COG_latent, DATE_cog = EXAMDATE), .(PTID, EXAMDATE)] |>
  unique()

#DT_c        <-
#DT_c        <- age.dt[cog.dt, on = .(PTID, VISCODE), .(PTID, AGE, COG_latent)
                      #][DT[, .(EXAMDATE, SEX, APOE4, HCv, HVR,
                               #WMHvol, AMYLOID, TAU, PATHOLOGY_lat, PET_lat,
                               #MRI_lat), .(PTID, AGE)], on = .(PTID, AGE)] |>
#na.omit() |>
#setkey(PTID)

#setcolorder(DT_c, c(1:2,4:6))
#DT_c[, COG_orig := COG_latent]
#DT_c[, COG_latent := COG_latent + abs(min(COG_latent))]
#DT_c        <- DTclean(DT_c, ordervars = "APOE4", scalevars = names(DT_c)[6:14])

#hcv   <- y
#y     <- "COG_latent_scl"
#names(mods) <- hcv
#base_mod    <- paste("COG_latent_scl ~ TIME + %s",
                     #"+ SEX + AGE.bl + APOE4 + (1|PTID)")
#base_mod2   <- sub("\\+", "*", base_mod)
#f2 <- c(sprintf(base_mod, hcv[1]),
        #sprintf(base_mod, hcv[2]),
        #sprintf(base_mod, "AMYLOID_scl"),
        #sprintf(base_mod, "TAU_scl"),
        #sprintf(base_mod, "PET_lat_scl"),
        #sprintf(base_mod, "PATHOLOGY_lat_scl"),
        #sprintf(base_mod, "MRI_lat_scl"),
        #sprintf(base_mod, "MRI_lat_scl + PET_lat_scl"),
        #sprintf(base_mod, "MRI_lat_scl + TAU_scl"))
#mods <- vector("list", length(f2))
#for (j in seq_along(f2)) mods[j] <- f2[j]

### Fits
#fpath <- here("data/rds/adni_mlm-cog-pathology.rds")
### Progress bar
#pb <- progress_bar$new(format = "MLMs | :what [:bar] :current/:total",
                       #total = length(y) * length(f2),
                       #clear = FALSE, width = 75)
#if (!file.exists(fpath) | ReFitModels) {
##if (TRUE) {
  #fits2  <- vector("list", length(mods))
  #for (j in seq_along(mods)) {
    ##pb$tick(tokens = list(what = sprintf("%s — model %i", y[i], j)))
    #fits2[[j]] <- lmer(as.formula(mods[[j]]), DT_c, REML = F,
                       #control = lmerControl(optimizer = "bobyqa"))
  #}
  #saveRDS(fits2, fpath)
#} else {
  #fits2         <- readRDS(fpath)
#}

### Observed variables
#obs_labels  <- c("Hippocampal Volume", "HVR",
                 #"Amyloid (PET)", "Tau (PET)")
#lat_labels  <- c("F: TAU + PET", "F: HVR + WMH", "F: WMH + TAU + PET")
#cov_labels  <- c("Intercept", "Time (Years)", "Sex (Male)",
                 #"Age (at baseline)", "APOE4 (Lin.)", "APOE4 (Quad.)")

#stargazer(fits2[1:4], intercept.bottom = F, type = "html",
          #out = here("data/derivatives/adni_mlm_cog-pathology1.html"),
          #title = "Mixed effect models (random intercept; standardized)",
          #dep.var.labels = "Cognitive Decline (Latent variable)",
          #covariate.labels = c(cov_labels[1:2], obs_labels, cov_labels[3:6]))

#stargazer(fits2[5:7], intercept.bottom = F, type = "html",
          #out = here("data/derivatives/adni_mlm_cog-pathology2.html"),
          #title = "Mixed effect models (random intercept; standardized)",
          #dep.var.labels = "Cognitive Decline (Latent variable)",
          #covariate.labels = c(cov_labels[1:2], lat_labels, cov_labels[3:6]))
