#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
library(stringr)
library(progress)
library(lavaan)
library(lavaanPlot)
library(ggplot2)
library(semTable)

## Refit models
refit_mods  <- FALSE

## Print plots ### Needs to be done outside renv
print_plots <- FALSE

## INPUT
# Load baseline data
fpath       <- here("data/rds/triad.rds")
if (file.exists(fpath)) {
  triad.dt  <- read_rds(fpath)
} else {
  here("code/demographics.R") |> source()
}
rm(fpath)

# Load Cerebra data
# Biomarkers
fpath       <- here("data/pet_biomarkers_cerebra.csv")
if (file.exists(fpath)) {
  cerebra.dt    <- fread(fpath)
} else {
  here("code/parse_csv_data.R") |> source()
}

# ROIs
fpaths      <- here("data/rds",
                    sprintf("cerebra_rois_%s_moca.rds",
                            c("amy", "tau", "amy_tau")))
if (any(!file.exists(fpaths))) {
  here("code/feature_selection.R") |> source()
} else {
  rois_amy      <- read_rds(fpaths[1])
  rois_tau      <- read_rds(fpaths[2])
  rois_amy_tau  <- read_rds(fpaths[3])
}
rm(fpaths)

## Data cleaning
# Convert Sex to dummy variable
triad.dt[, SEX_n := as.numeric(SEX) - 1]

# Remove youth and dementias
triad.dt    <- triad.dt[!DX_clean %in% c("Young", "Other", "AD")]

# 1 - HVR (average for both sides)
#triad.dt[, `:=`(HVR_lr = 1 - HVR_l, HVR_rr = 1 - HVR_r)]
triad.dt[, `:=`(HVR_mean_inv = 1 - (HVR_l + HVR_r) / 2)]

triad.dt    <- triad.dt[, .(PTID, VISIT, DX, SEX_n, AGE_scan, EDUC,#APOE_n,
                            HVR_mean_inv, MOCA_score)]

# Merge triad and cerebra data
# Cerebra dictionary
dict_roi    <- unique(cerebra.dt[, .(LABEL_id, LABEL_name, SIDE)])

# Selected features for both Amyloid and Tau
amy.dt      <- cerebra.dt[LABEL_id %in% rois_amy_tau$LABEL_id
                          ][!is.na(AMYLOID_norm),
                          .(PTID, VISIT, LABEL_id, AMYLOID = AMYLOID_norm)]

tau.dt      <- cerebra.dt[LABEL_id %in% rois_amy_tau$LABEL_id
                          ][!is.na(AMYLOID_norm),
                          .(PTID, VISIT, LABEL_id, TAU = TAU_norm)]

triad.dt    <- tau.dt[amy.dt, on = .(PTID, VISIT, LABEL_id)
                      ][triad.dt, on = .(PTID, VISIT)
                      ][!is.na(MOCA_score) & !is.na(AMYLOID) & !is.na(TAU)]

## Model definitions
hvr.mod     <-
  str_glue("
# Regressions
AMYLOID ~ SEX_n + AGE_scan
TAU ~ a * AMYLOID + SEX_n + AGE_scan
HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan
# HVR mediation
# Direct effect
dAMY  := b
# Indirect effect
iTAU  := a * c
# Total effect
Total := dAMY + iTAU
# Proportion analysis
propAMY := dAMY / Total
propTAU := iTAU / Total
")

moca.mod    <-
  str_glue("
# Regressions
AMYLOID ~ SEX_n + AGE_scan
TAU ~ a * AMYLOID + SEX_n + AGE_scan
HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan
MOCA_score ~ d * AMYLOID + e * TAU + f * HVR_mean_inv + SEX_n + AGE_scan + EDUC
# Direct effect
dAMY := d
# Indirect effects
iTAU := a * e
iHVR := (b + (a * c)) * f
# Total effects
iTotal    := iTAU + iHVR
Total := iTotal + dAMY
# Proportion analysis
propAMY := dAMY / Total
propTAU := iTAU / Total
propHVR := iHVR / Total
propInd := iTotal / Total
")


## Fit models
fnames          <- here("data/rds",
                        sprintf("mediation_%s_bs.rds", c("hvr", "moca")))

if (any(!file.exists(fnames)) | refit_mods) {
  mod_hvr.fits  <- mod_moca.fits <- vector("list", rois_amy_tau[, .N])
  mods.lst      <- c(mod_hvr.fits, mod_moca.fits)

  pb <- progress_bar$new(format = "Models | :what [:bar] :current/:total",
                         total = length(mods.lst) * rois_amy_tau[, .N],
                         clear = FALSE, width = 75)

  for (i in seq_along(mods.lst)) {
    names(mods.lst[i])  <- rois_amy_tau[, paste(LABEL_name, SIDE, sep = "_")]

    for (j in seq_along(rois_amy_tau$LABEL_id)) {
      pb$tick(tokens = list(what = sprintf("%s : %s",
                                           c("HVR", "MoCA")[i],
                                           rois_amy_tau[j, LABEL_id])))

      mods.lst[i][[j]]  <-
        sem(c(hvr.mod, moca.mod)[i],
            data = triad.dt[LABEL_id == rois_amy_tau[i, LABEL_id]],
            estimator = "ML",
            se = "bootstrap",
            bootstrap = 10000)
    }

    write_rds(mods.lst[i], fnames[i])
  }
} else {
  mod_hvr.fits          <- read_rds(fnames[1])
  mod_moca.fits         <- read_rds(fnames[2])
}
rm(fnames)

## Extract Fit measures and paramater estimates
# Names for data cleaning
msrs_names      <- mod_hvr.fits[[1]] |> fitMeasures() |> names()
rois_names      <- rois_amy_tau[, paste(LABEL_name, SIDE, sep = "_")]

# HVR model
mod_hvr.msrs    <- mod_hvr.fits |>
  lapply(fitMeasures) |>
  lapply(as.data.table) |>
  lapply(transpose) |>
  rbindlist()

setnames(mod_hvr.msrs, msrs_names)
mod_hvr.msrs[, ROI := rois_names]

mod_hvr.est     <- mod_hvr.fits |>
  lapply(standardizedSolution) |>
  lapply(as.data.table)

mod_hvr.est     <- rois_names |>
  lapply(function (name) mod_hvr.est[[name]][, ROI := name]) |>
  rbindlist()

mod_hvr.est     <- mod_hvr.est[op != "~~"]

mod_moca.msrs   <- mod_moca.fits |>
  lapply(fitMeasures) |>
  lapply(as.data.table) |>
  lapply(transpose) |>
  rbindlist()

setnames(mod_moca.msrs, msrs_names)
mod_hvr.msrs[, ROI := rois_names]

mod_moca.est    <- mod_moca.fits |>
  lapply(standardizedSolution) |>
  lapply(as.data.table)

mod_moca.est    <- rois_names |>
  lapply(function (name) mod_moca.est[[name]][, ROI := name]) |>
  rbindlist()

mod_moca.est    <- mod_moca.est[op != "~~"]

## Plot standardized estimates
# HVR model
mod_hvr.est[lhs == "dAMY", label := "Direct: Amyloid"]
mod_hvr.est[lhs == "iTAU", label := "Indirect: Tau"]

mod_hvr.est[lhs %in% c("dAMY", "iTAU")] |>
  ggplot(aes(x = ROI, y = est.std)) +
  theme_classic(base_size = 12) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper),
                width = 0.2, position = position_dodge(0.5)) +
  geom_point(shape = 21, fill = "white", size = 1.5) +
  labs(x = "CerebrA ROIs", y = "Standardized estimates",
       title = "Direct & indirect effects on HC atrophy",
       caption = "Bootstrap CIs: 10000 resamples") +
  coord_flip() +
  facet_wrap(vars(label))

here("plots/mediation_hvr_std-estimates.png") |>
  ggsave(width = 7, height = 5, units = "in", dpi = 600)

# MoCA model
mod_moca.est[lhs == "dAMY", label := "Direct: Amyloid"]
mod_moca.est[lhs == "iTAU", label := "Indirect: Tau"]
mod_moca.est[lhs == "iHVR", label := "Indirect: HVR"]

mod_moca.est[lhs %in% c("dAMY", "iTAU")] |>
  ggplot(aes(x = ROI, y = est.std)) +
  theme_classic(base_size = 12) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper),
                width = 0.2, position = position_dodge(0.5)) +
  geom_point(shape = 21, fill = "white", size = 1.5) +
  labs(x = "CerebrA ROIs", y = "Standardized estimates",
       title = "Direct & indirect effects on MoCA scores",
       caption = "Bootstrap CIs: 10000 resamples") +
  coord_flip() +
  facet_wrap(vars(label))

here("plots/mediation_moca_std-estimates2.png") |>
  ggsave(width = 7, height = 5, units = "in", dpi = 600)


## RMSEA = 0; CFI & TLI = 1 #
#fname <- here("data/rds/mediation_selected_moca.rds")
#if (!file.exists(fname) | refit_mods) {
  #mod_sel_moca.fit  <- sem(moca.mod, data = triad_s.dt, #group = "DX",
                                 #estimator = "ML")
                     ##se = "bootstrap", bootstrap = 10000)
  #write_rds(mod_sel_moca.fit, fname)
#} else {
  #mod_sel_moca.fit  <- read_rds(fname)
#}
#rm(fname)

## RMSEA > .18; low CFI/TLI
#fname <- here("data/rds/mediation_latent_hvr.rds")
#if (!file.exists(fname) | refit_mods) {
  #mod_lat_hvr.fit <- sem(hvr_kl.mod, data = triad_k.dt, #group = "DX",
                                #estimator = "ML")
                     ##se = "bootstrap", bootstrap = 10000)
  #write_rds(mod_lat_hvr.fit, fname)
#} else {
  #mod_lat_hvr.fit <- read_rds(fname)
#}
#rm(fname)

## RMSEA > .16; low CFI/TLI
#fname <- here("data/rds/mediation_latent_moca.rds")
#if (!file.exists(fname) | refit_mods) {
  #mod_lat_moca.fit <- sem(moca_kl.mod, data = triad_k.dt, #group = "DX",
                                 #estimator = "ML")
                     ##se = "bootstrap", bootstrap = 10000)
  #write_rds(mod_lat_moca.fit, fname)
#} else {
  #mod_lat_moca.fit <- read_rds(fname)
#}
#rm(fname)

## Unusable: RMSEA > .5; CFI < .3; TLI < 0
#fname <- here("data/rds/mediation_clusters_hvr.rds")
#if (!file.exists(fname) | refit_mods) {
  #mod_clust_hvr.fit <- sem(hvr_kd.mod, data = triad_k.dt, #group = "DX",
                                #estimator = "ML")
                     ##se = "bootstrap", bootstrap = 10000)
  #write_rds(mod_clust_hvr.fit, fname)
#} else {
  #mod_clust_hvr.fit <- read_rds(fname)
#}
#rm(fname)

## Unusable: RMSEA ~ .5; CFI < .32; TLI < 0
#fname <- here("data/rds/mediation_clusters_moca.rds")
#if (!file.exists(fname) | refit_mods) {
  #mod_clust_moca.fit <- sem(moca_kd.mod, data = triad_k.dt, #group = "DX",
                                 #estimator = "ML")
                     ##se = "bootstrap", bootstrap = 10000)
  #write_rds(mod_clust_moca.fit, fname)
#} else {
  #mod_clust_moca.fit <- read_rds(fname)
#}
#rm(fname)

#### Simplest model with ALL ROIs
### HVR mediation
##fname <- here("data/rds/mediation_simplest_hvr.rds")
##if (!file.exists(fname) | refit_mods) {
  ##model_simplest_hvr.fit <- sem(hvr.mod,
                                ##data = triad_all.dt,
                                ###group = "DX",
                                ##estimator = "ML")
                     ###se = "bootstrap", bootstrap = 10000)
  ##write_rds(model_simplest_hvr.fit, fname)
  ##semTable(model_simplest_hvr.fit, print.results = FALSE,
           ##paramSets  = c("slopes", "constructed", "fits"),
           ##columns    = c("estsestars", "z"),
           ##fits       = c("npar", "chisq", "pvalue", "cfi", "tli",
                          ##"aic", "bic2", "rmsea"),
           ##type       = "html",
           ##file       = here("data/derivatives/mediation_simplest_hvr.html"))
##} else {
  ##model_simplest_hvr.fit <- read_rds(fname)
##}
##rm(fname)

### MoCA mediation
##fname <- here("data/rds/mediation_rois_simplest_moca.rds")
##if (!file.exists(fname) | refit_mods) {
  ##model_simplest_moca.fit <- sem(moca.mod,
                                 ##data = triad_all.dt,
                                 ###group = "DX",
                                 ##estimator = "ML")
                     ###se = "bootstrap", bootstrap = 10000)
  ##write_rds(model_simplest_moca.fit, fname)
  ##semTable(model_simplest_moca.fit, print.results = FALSE,
           ##paramSets  = c("slopes", "constructed", "fits"),
           ##columns    = c("estsestars", "z"),
           ##fits       = c("npar", "chisq", "pvalue", "cfi", "tli",
                          ##"aic", "bic2", "rmsea"),
           ##type       = "html",
           ##file       = here("data/derivatives/mediation_simplest_moca.html"))
##} else {
  ##model_simplest_moca.fit <- read_rds(fname)
##}
##rm(fname)

#### Simple model with selected ROIs
### HVR mediation
##fname <- here("data/rds/mediation_simple_hvr.rds")
##if (!file.exists(fname) | refit_mods) {
  ##model_simple_hvr.fit <- sem(hvr.mod,
                               ##data = triad_sel.dt,
                               ###group = "DX",
                               ##estimator = "ML")
                     ###se = "bootstrap", bootstrap = 10000)
  ##write_rds(model_simple_hvr.fit, fname)
  ##semTable(model_simple_hvr.fit, print.results = FALSE,
           ##paramSets  = c("slopes", "constructed", "fits"),
           ##columns    = c("estsestars", "z"),
           ##fits       = c("npar", "chisq", "pvalue", "cfi", "tli",
                          ##"aic", "bic2", "rmsea"),
           ##type       = "html",
           ##file       = here("data/derivatives/mediation_simple_hvr.html"))
##} else {
  ##model_simple_hvr.fit <- read_rds(fname)
##}
##rm(fname)

### MoCA mediation
##fname <- here("data/rds/mediation_rois_simple_moca.rds")
##if (!file.exists(fname) | refit_mods) {
  ##model_simple_moca.fit <- sem(moca.mod,
                               ##data = triad_sel.dt,
                               ###group = "DX",
                               ##estimator = "ML")
                     ###se = "bootstrap", bootstrap = 10000)
  ##write_rds(model_simple_moca.fit, fname)
  ##semTable(model_simple_moca.fit, print.results = FALSE,
           ##paramSets  = c("slopes", "constructed", "fits"),
           ##columns    = c("estsestars", "z"),
           ##fits       = c("npar", "chisq", "pvalue", "cfi", "tli",
                          ##"aic", "bic2", "rmsea"),
           ##type       = "html",
           ##file       = here("data/derivatives/mediation_simple_moca.html"))
##} else {
  ##model_simple_moca.fit <- read_rds(fname)
##}
##rm(fname)

#### Labels
##labels_cov  <- c(AGE_scan = "Age", SEX_n = "Sex", APOE_n = "APOE4")
##labels_hcv  <- c(HCv_l = "Left", HCv_r = "Right")
##labels_hvr  <- c(HVR = "HC-atrophy",
                 ##HVR_lr = "1-HVR (Left)",
                 ##HVR_rr = "1-HVR (Right)")
##labels_hvr2 <- c(HVR_mean_inv = "HC-atrophy")
##labels_moca <- c(MOCA_score = "MoCA")

#### Models with latent variables ###
##amy_ftrs  <- paste("AMY", 1:5, sep = "_")
##tau_ftrs  <- paste("TAU", 1:5, sep = "_")

##latent.mod <-
  ##str_glue("
           ### Latent variables
           ##AMYLOID =~ {paste(amy_ftrs, collapse = ' + ')}
           ##TAU =~ {paste(tau_ftrs, collapse = ' + ')}
           ### Regressions
           ##AMYLOID ~ SEX_n + AGE_scan
           ##TAU ~ a * AMYLOID + SEX_n + AGE_scan
           ##HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan
           ##MOCA_score ~ d * AMYLOID + e * TAU + f * HVR_mean_inv + SEX_n + AGE_scan + EDUC
           ### HVR mediation
           ### Direct effect
           ##dAMY_HVR  := b
           ### Indirect effect
           ##iTAU_HVR  := a * c
           ### Total effect
           ##Total_HVR := dAMY_HVR + iTAU_HVR
           ### MOCA mediation
           ### Direct effect
           ##dAMY_MOCA := d
           ### Indirect effects
           ##iTAU_MOCA := a * e
           ##iHVR_MOCA := Total_HVR * f
           ### Total effects
           ##iTotal    := iTAU_MOCA + iHVR_MOCA
           ##Total_MOCA := iTotal + dAMY_MOCA
           ### Proportion analysis
           ##propAMY_HVR := dAMY_HVR / Total_HVR
           ##propTAU_HVR := iTAU_HVR / Total_HVR
           ##propAMY_MOCA := dAMY_MOCA / Total_MOCA
           ##propTAU_MOCA := iTAU_MOCA / Total_MOCA
           ##propHVR_MOCA := iHVR_MOCA / Total_MOCA
           ##propIndirect_MOCA := iTotal / Total_MOCA
           ##")

##fname <- here("data/rds/mediation_rois_latent.rds")
##if (!file.exists(fname) | refit_mods) {
  ##model_rois_latent.fit <- sem(latent.mod,
                               ##data = triad.dt[order(DX)],
                               ###group = "DX",
                               ##estimator = "ML")
                     ###se = "bootstrap", bootstrap = 10000)
  ##write_rds(model_rois_latent.fit, fname)
##} else {
  ##model_rois_latent.fit <- read_rds(fname)
##}
##rm(fname)

##### Model with all mediators ###
### Text parsing
##amy_covars    <- paste(amy_ftrs, collapse = " + ")
##tau_covars    <- paste(tau_ftrs, collapse = " + ")

##amy_regress   <- sprintf("%s ~ SEX_n + AGE_scan", amy_ftrs)
##tau_regress   <- sprintf("%s ~ %s + SEX_n + AGE_scan",
                         ##tau_ftrs, amy_covars)

### Model definition
##detailed.mod <-
  ##str_glue("
#### Regressions ##
### Amyloid
##{paste(amy_regress, collapse = '\n')}
### Tau
##{paste(tau_regress, collapse = '\n')}
##HVR_mean_inv ~ {amy_covars} + {tau_covars} + SEX_n + AGE_scan
##MOCA_score ~ {amy_covars} + {tau_covars} + HVR_mean_inv + SEX_n + AGE_scan + EDUC
##")

##fname <- here("data/rds/mediation_rois_detailed.rds")
##if (!file.exists(fname) | refit_mods) {
  ##model_rois_detailed.fit <- sem(detailed.mod,
                                 ##data = triad.dt[order(DX)],
                                 ###group = "DX",
                                 ##estimator = "ML")
                     ###se = "bootstrap", bootstrap = 10000)
  ##write_rds(model_rois_detailed.fit, fname)
##} else {
  ##model_rois_detailed.fit <- read_rds(fname)
##}
##rm(fname)

### Estimates data.table
##estimates.dt  <- parameterEstimates(model_rois_detailed.fit,
                                    ##standardized = TRUE) |>
           ##as.data.table()

### Remove variances data
##estimates.dt  <- estimates.dt[op != "~~"]

### Adjust for multiple comparisons
##estimates.dt[, pvalue_adj := p.adjust(pvalue, method = "bonferroni")]

##estimates.dt[, lhs_roi_id := as.numeric(str_extract(lhs, "\\d{3}"))]
##estimates.dt[, rhs_roi_id := as.numeric(str_extract(rhs, "\\d{3}"))]
