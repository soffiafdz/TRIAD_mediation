#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
library(stringr)
library(lavaan)
library(lavaanPlot)

## Refit models
refit_models <- TRUE

## Print plots ### Needs to be done outside renv
print_plots <- FALSE
#print_plots <- TRUE

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
  cerebra.dt  <- fread(fpath)
} else {
  here("code/parse_csv_data.R") |> source()
}

# ROIs
fpaths      <- here("data/rds",
                    sprintf("cerebra_rois_%s_moca.rds", c("amy", "tau")))
if (any(!file.exists(fpaths))) {
  here("code/feature_selection.R") |> source()
} else {
  rois_amy  <- read_rds(fpaths[1])
  rois_tau  <- read_rds(fpaths[2])
}
rm(fpaths)

## Data cleaning
# Convert Sex to dummy variable
triad.dt[, SEX_n := as.numeric(SEX) - 1]

# Remove youth and other dementias
triad.dt    <- triad.dt[!DX_clean %in% c("Young", "Other")]
triad.dt[, DX := factor(DX, levels = c("CN", "MCI", "AD"))]

# 1 - HVR (average for both sides)
#triad.dt[, `:=`(HVR_lr = 1 - HVR_l, HVR_rr = 1 - HVR_r)]
triad.dt[, `:=`(HVR_mean_inv = 1 - (HVR_l + HVR_r) / 2)]

triad.dt    <- triad.dt[, .(PTID, VISIT, TAU_braak_group,
                            SEX_n, AGE_scan, EDUC, APOE_n,
                            HVR_mean_inv, MOCA_score)]

# Merge triad and cerebra data
# Cerebra dictionary
dict_roi    <- unique(cerebra.dt[, .(LABEL_id, LABEL_name, SIDE)])

# Amyloid
amy_roi.dt  <- cerebra.dt[!is.na(AMYLOID),
                          .(PTID, VISIT, AMYLOID,
                            ROI = sprintf("AMY_%03i", LABEL_id))] |>
  dcast(... ~ ROI, value.var = "AMYLOID")

amy_roi.dt  <- amy_roi.dt[triad.dt, on = .(PTID, VISIT)]

# Tau
tau_roi.dt  <- cerebra.dt[!is.na(TAU),
                          .(PTID, VISIT, TAU,
                            ROI = sprintf("TAU_%03i", LABEL_id))] |>
  dcast(... ~ ROI, value.var = "TAU")

triad.dt  <- tau_roi.dt[amy_roi.dt, on = .(PTID, VISIT)]
triad.dt  <- triad.dt[!is.na(AMY_001) &
                      !is.na(TAU_001) &
                      !is.na(TAU_braak_group)]
rm(amy_roi.dt, tau_roi.dt)

## Labels
labels_cov  <- c(AGE_scan = "Age", SEX_n = "Sex", APOE_n = "APOE4")
labels_hcv  <- c(HCv_l = "Left", HCv_r = "Right")
labels_hvr  <- c(HVR = "HC-atrophy",
                 HVR_lr = "1-HVR (Left)",
                 HVR_rr = "1-HVR (Right)")
labels_tau  <- c(TAU_braak1 = "Braak1",
                 TAU_braak2 = "Braak2",
                 TAU_braak3 = "Braak3",
                 TAU_braak4 = "Braak4",
                 TAU_braak5 = "Braak5",
                 TAU_braak6 = "Braak6")
labels_hvr2 <- c(HVR_mean_inv = "HC-atrophy")
labels_moca <- c(MOCA_score = "MoCA")

### Models with latent variables ###
# Lack of Tau aggregation (Braak Stage 0)
# Confirmed important features from Boruta algorithm
amy_features  <- rois_amy[[1]][decision == "Confirmed", id]
tau_features  <- rois_tau[[1]][decision == "Confirmed", id]

latent1.mod <-
  str_glue("
           # Latent variables
           AMYLOID =~ {paste(amy_features, collapse = ' + ')}
           TAU =~ {paste(tau_features, collapse = ' + ')}
           # Regressions
           AMYLOID ~ SEX_n + AGE_scan
           TAU ~ a * AMYLOID + SEX_n + AGE_scan
           HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan
           MOCA_score ~ d * AMYLOID + e * TAU + f * HVR_mean_inv + SEX_n + AGE_scan + EDUC
           # Total effect
           Total := d + (a * e) + (b * f)
           # Indirect effects
           ieAMY := d
           propAMY := ieAMY / Total
           ieTAU := a * e
           propTAU := ieTAU / Total
           ieHVR := b * f
           propHVR := ieHVR / Total
           ieTotal := (a * e) + (b * f)
           propTotal := ieTotal/ Total
           ")

fname <- here("data/rds/mediation_rois_latent1.rds")
if (!file.exists(fname) | refit_models) {
  model_rois_latent1.fit <- sem(latent1.mod,
                               data = triad.dt[TAU_braak_group == "0"],
                               estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(model_rois_latent1.fit, fname)
} else {
  model_rois_latent1.fit <- read_rds(fname)
}
rm(fname)

# Initial stages of Tau aggregation (Braak Stage 1 & 2)
# Confirmed important features from Boruta algorithm
amy_features  <- rois_amy[[2]][decision == "Confirmed", id]
tau_features  <- rois_tau[[2]][decision == "Confirmed", id]

latent2.mod <-
  str_glue("
           # Latent variables
           AMYLOID =~ {paste(amy_features, collapse = ' + ')}
           TAU =~ {paste(tau_features, collapse = ' + ')}
           # Regressions
           AMYLOID ~ SEX_n + AGE_scan
           TAU ~ a * AMYLOID + SEX_n + AGE_scan
           HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan
           MOCA_score ~ d * AMYLOID + e * TAU + f * HVR_mean_inv + SEX_n + AGE_scan + EDUC
           # Total effect
           Total := d + (a * e) + (b * f)
           # Indirect effects
           ieAMY := d
           propAMY := ieAMY / Total
           ieTAU := a * e
           propTAU := ieTAU / Total
           ieHVR := b * f
           propHVR := ieHVR / Total
           ieTotal := (a * e) + (b * f)
           propTotal := ieTotal/ Total
           ")

fname <- here("data/rds/mediation_rois_latent2.rds")
if (!file.exists(fname) | refit_models) {
  model_rois_latent2.fit <- sem(latent2.mod,
                               data = triad.dt[TAU_braak_group == "1 & 2"],
                               estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(model_rois_latent2.fit, fname)
} else {
  model_rois_latent2.fit <- read_rds(fname)
}
rm(fname)

# Moderate stages of Tau aggregation (Braak Stage 3 & 4)
# Confirmed important features from Boruta algorithm
amy_features  <- rois_amy[[3]][decision == "Confirmed", id]
tau_features  <- rois_tau[[3]][decision == "Confirmed", id]

latent3.mod <-
  str_glue("
           # Latent variables
           AMYLOID =~ {paste(amy_features, collapse = ' + ')}
           TAU =~ {paste(tau_features, collapse = ' + ')}
           # Regressions
           AMYLOID ~ SEX_n + AGE_scan
           TAU ~ a * AMYLOID + SEX_n + AGE_scan
           HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan
           MOCA_score ~ d * AMYLOID + e * TAU + f * HVR_mean_inv + SEX_n + AGE_scan + EDUC
           # Total effect
           Total := d + (a * e) + (b * f)
           # Indirect effects
           ieAMY := d
           propAMY := ieAMY / Total
           ieTAU := a * e
           propTAU := ieTAU / Total
           ieHVR := b * f
           propHVR := ieHVR / Total
           ieTotal := (a * e) + (b * f)
           propTotal := ieTotal/ Total
           ")

fname <- here("data/rds/mediation_rois_latent3.rds")
if (!file.exists(fname) | refit_models) {
  model_rois_latent3.fit <- sem(latent3.mod,
                               data = triad.dt[TAU_braak_group == "3 & 4"],
                               estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(model_rois_latent3.fit, fname)
} else {
  model_rois_latent3.fit <- read_rds(fname)
}
rm(fname)

# Severe stages of Tau aggregation (Braak Stage 5 & 6)
# Confirmed important features from Boruta algorithm
amy_features  <- rois_amy[[4]][decision == "Confirmed", id]
tau_features  <- rois_tau[[4]][decision == "Confirmed", id]

latent4.mod <-
  str_glue("
           # Latent variables
           AMYLOID =~ {paste(amy_features, collapse = ' + ')}
           TAU =~ {paste(tau_features, collapse = ' + ')}
           # Regressions
           AMYLOID ~ SEX_n + AGE_scan
           TAU ~ a * AMYLOID + SEX_n + AGE_scan
           HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan
           MOCA_score ~ d * AMYLOID + e * TAU + f * HVR_mean_inv + SEX_n + AGE_scan + EDUC
           # Total effect
           Total := d + (a * e) + (b * f)
           # Indirect effects
           ieAMY := d
           propAMY := ieAMY / Total
           ieTAU := a * e
           propTAU := ieTAU / Total
           ieHVR := b * f
           propHVR := ieHVR / Total
           ieTotal := (a * e) + (b * f)
           propTotal := ieTotal/ Total
           ")

fname <- here("data/rds/mediation_rois_latent4.rds")
if (!file.exists(fname) | refit_models) {
  model_rois_latent4.fit <- sem(latent4.mod,
                               data = triad.dt[TAU_braak_group == "5 & 6"],
                               estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(model_rois_latent4.fit, fname)
} else {
  model_rois_latent4.fit <- read_rds(fname)
}
rm(fname)

### Model with all mediators ###
## Lack of Tau aggregation (Braak Stage 0)
amy_features  <- rois_amy[[4]][decision == "Confirmed", id]
tau_features  <- rois_tau[[4]][decision == "Confirmed", id]

# Text parsing
amy_covars    <- paste(amy_features, collapse = " + ")
tau_covars    <- paste(tau_features, collapse = " + ")

amy_regress   <- sprintf("%s ~ SEX_n + AGE_scan", amy_features)
tau_regress   <- sprintf("%s ~ %s + SEX_n + AGE_scan",
                         tau_features, amy_covars)

# Model definition
detailed1.mod <-
  str_glue("
## Regressions ##
# Amyloid
{paste(amy_regress, collapse = '\n')}
# Tau
{paste(tau_regress, collapse = '\n')}
HVR_mean_inv ~ {amy_covars} + {tau_covars} + SEX_n + AGE_scan
MOCA_score ~ {amy_covars} + {tau_covars} + HVR_mean_inv + SEX_n + AGE_scan + EDUC
")

fname <- here("data/rds/mediation_rois_detailed1.rds")
if (!file.exists(fname) | refit_models) {
  model_rois_detailed1.fit <- sem(detailed1.mod,
                               data = triad.dt,
                               estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(model_rois_detailed1.fit, fname)
} else {
  model_rois_detailed1.fit <- read_rds(fname)
}
rm(fname)
