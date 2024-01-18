#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
library(stringr)
library(lavaan)
library(lavaanPlot)

## Refit models
#refit_models <- TRUE
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
triad.dt    <- triad.dt[!DX_clean %in% c("Young", "Other", "AD")]
#triad.dt    <- triad.dt[!DX_clean %in% c("Young", "Other")]

# 1 - HVR (average for both sides)
#triad.dt[, `:=`(HVR_lr = 1 - HVR_l, HVR_rr = 1 - HVR_r)]
triad.dt[, `:=`(HVR_mean_inv = 1 - (HVR_l + HVR_r) / 2)]

triad.dt    <- triad.dt[, .(PTID, VISIT, DX, SEX_n, AGE_scan, EDUC,#APOE_n,
                            HVR_mean_inv, MOCA_score)]

# Merge triad and cerebra data
# Cerebra dictionary
dict_roi    <- unique(cerebra.dt[, .(LABEL_id, LABEL_name, SIDE)])

# Amyloid
amy.dt      <- rois_amy.dt[!is.na(cluster),
                           .(LABEL_id = as.numeric(str_extract(id, "\\d{3}")),
                             CLUSTER  = sprintf("AMY_%i", cluster))
                           ][cerebra.dt, on = "LABEL_id"
                           ][!is.na(AMYLOID) & !is.na(CLUSTER),
                           weighted.mean(AMYLOID, VOL),
                           .(PTID, VISIT, CLUSTER)] |>
  dcast(... ~ CLUSTER, value.var = "V1")

amy.dt  <- amy.dt[triad.dt, on = .(PTID, VISIT)]

# Tau
tau.dt      <- rois_tau.dt[!is.na(cluster),
                           .(LABEL_id = as.numeric(str_extract(id, "\\d{3}")),
                             CLUSTER  = sprintf("TAU_%i", cluster))
                           ][cerebra.dt, on = "LABEL_id"
                           ][!is.na(TAU) & !is.na(CLUSTER),
                           weighted.mean(TAU, VOL),
                           .(PTID, VISIT, CLUSTER)] |>
  dcast(... ~ CLUSTER, value.var = "V1")

triad.dt  <- tau.dt[amy.dt, on = .(PTID, VISIT)]
triad.dt  <- triad.dt[!is.na(AMY_1) & !is.na(TAU_1) & !is.na(MOCA_score)]
rm(amy.dt, tau.dt)

### Labels
#labels_cov  <- c(AGE_scan = "Age", SEX_n = "Sex", APOE_n = "APOE4")
#labels_hcv  <- c(HCv_l = "Left", HCv_r = "Right")
#labels_hvr  <- c(HVR = "HC-atrophy",
                 #HVR_lr = "1-HVR (Left)",
                 #HVR_rr = "1-HVR (Right)")
#labels_hvr2 <- c(HVR_mean_inv = "HC-atrophy")
#labels_moca <- c(MOCA_score = "MoCA")

### Models with latent variables ###
amy_features  <- paste("AMY", 1:5, sep = "_")
tau_features  <- paste("TAU", 1:5, sep = "_")

latent.mod <-
  str_glue("
           # Latent variables
           AMYLOID =~ {paste(amy_features, collapse = ' + ')}
           TAU =~ {paste(tau_features, collapse = ' + ')}
           # Regressions
           AMYLOID ~ SEX_n + AGE_scan
           TAU ~ a * AMYLOID + SEX_n + AGE_scan
           HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan
           MOCA_score ~ d * AMYLOID + e * TAU + f * HVR_mean_inv + SEX_n + AGE_scan + EDUC
           # HVR mediation
           # Direct effect
           dAMY_HVR  := b
           # Indirect effect
           iTAU_HVR  := a * c
           # Total effect
           Total_HVR := dAMY_HVR + iTAU_HVR
           # MOCA mediation
           # Direct effect
           dAMY_MOCA := d
           # Indirect effects
           iTAU_MOCA := a * e
           iHVR_MOCA := Total_HVR * f
           # Total effects
           iTotal    := iTAU_MOCA + iHVR_MOCA
           Total_MOCA := iTotal + dAMY_MOCA
           # Proportion analysis
           propAMY_HVR := dAMY_HVR / Total_HVR
           propTAU_HVR := iTAU_HVR / Total_HVR
           propAMY_MOCA := dAMY_MOCA / Total_MOCA
           propTAU_MOCA := iTAU_MOCA / Total_MOCA
           propHVR_MOCA := iHVR_MOCA / Total_MOCA
           propIndirect_MOCA := iTotal / Total_MOCA
           ")

fname <- here("data/rds/mediation_rois_latent.rds")
if (!file.exists(fname) | refit_models) {
  model_rois_latent.fit <- sem(latent.mod,
                               data = triad.dt[order(DX)],
                               #group = "DX",
                               estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(model_rois_latent.fit, fname)
} else {
  model_rois_latent.fit <- read_rds(fname)
}
rm(fname)

### Model with all mediators ###
# Text parsing
amy_covars    <- paste(amy_features, collapse = " + ")
tau_covars    <- paste(tau_features, collapse = " + ")

amy_regress   <- sprintf("%s ~ SEX_n + AGE_scan", amy_features)
tau_regress   <- sprintf("%s ~ %s + SEX_n + AGE_scan",
                         tau_features, amy_covars)

# Model definition
detailed.mod <-
  str_glue("
## Regressions ##
# Amyloid
{paste(amy_regress, collapse = '\n')}
# Tau
{paste(tau_regress, collapse = '\n')}
HVR_mean_inv ~ {amy_covars} + {tau_covars} + SEX_n + AGE_scan
MOCA_score ~ {amy_covars} + {tau_covars} + HVR_mean_inv + SEX_n + AGE_scan + EDUC
")

fname <- here("data/rds/mediation_rois_detailed.rds")
if (!file.exists(fname) | refit_models) {
  model_rois_detailed.fit <- sem(detailed.mod,
                                 data = triad.dt[order(DX)],
                                 #group = "DX",
                                 estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(model_rois_detailed.fit, fname)
} else {
  model_rois_detailed.fit <- read_rds(fname)
}
rm(fname)

# Estimates data.table
estimates.dt  <- parameterEstimates(model_rois_detailed.fit,
                                    standardized = TRUE) |>
           as.data.table()

# Remove variances data
estimates.dt  <- estimates.dt[op != "~~"]

# Adjust for multiple comparisons
estimates.dt[, pvalue_adj := p.adjust(pvalue, method = "bonferroni")]

estimates.dt[, lhs_roi_id := as.numeric(str_extract(lhs, "\\d{3}"))]
estimates.dt[, rhs_roi_id := as.numeric(str_extract(rhs, "\\d{3}"))]
