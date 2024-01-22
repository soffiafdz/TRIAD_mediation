#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
library(stringr)
library(lavaan)
library(lavaanPlot)
library(semTable)

## Refit models
refit_mods  <- TRUE

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
                    sprintf("cerebra_rois_%s_moca.rds", c("amy", "tau")))
if (any(!file.exists(fpaths))) {
  here("code/feature_selection.R") |> source()
} else {
  rois_amy      <- read_rds(fpaths[1])
  rois_tau      <- read_rds(fpaths[2])
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

# Weighted average of all selected features
amy.dt      <-
  cerebra.dt[LABEL_id %in% rois_amy.dt[decision == "Confirmed",
                                       as.numeric(str_extract(id, "\\d{3}"))]
             ][!is.na(AMYLOID_norm),
               .(AMYLOID = weighted.mean(AMYLOID_norm, VOL)),
               .(PTID, VISIT)]

tau.dt      <-
  cerebra.dt[LABEL_id %in% rois_tau.dt[decision == "Confirmed",
                                       as.numeric(str_extract(id, "\\d{3}"))]
             ][!is.na(TAU_norm),
               .(TAU = weighted.mean(TAU_norm, VOL)),
               .(PTID, VISIT)]

triad_s.dt  <- tau.dt[amy.dt, on = .(PTID, VISIT)
                      ][triad.dt, on = .(PTID, VISIT)
                      ][!is.na(MOCA_score) & !is.na(AMYLOID) & !is.na(TAU)]
rm(amy.dt, tau.dt)

# Weighted average of selected features by cluster
K_amy       <- rois_amy.dt[!is.na(cluster), max(cluster)]
K_tau       <- rois_tau.dt[!is.na(cluster), max(cluster)]

amy.dt      <- rois_amy.dt[!is.na(cluster),
                           .(LABEL_id = as.numeric(str_extract(id, "\\d{3}")),
                             CLUSTER  = sprintf("AMY_%i", cluster))
                           ][cerebra.dt, on = "LABEL_id"
                           ][!is.na(AMYLOID_norm) & !is.na(CLUSTER),
                           weighted.mean(AMYLOID_norm, VOL),
                           .(PTID, VISIT, CLUSTER)] |>
  dcast(... ~ CLUSTER, value.var = "V1")

tau.dt      <- rois_tau.dt[!is.na(cluster),
                           .(LABEL_id = as.numeric(str_extract(id, "\\d{3}")),
                             CLUSTER  = sprintf("TAU_%i", cluster))
                           ][cerebra.dt, on = "LABEL_id"
                           ][!is.na(TAU_norm) & !is.na(CLUSTER),
                           weighted.mean(TAU_norm, VOL),
                           .(PTID, VISIT, CLUSTER)] |>
  dcast(... ~ CLUSTER, value.var = "V1")

triad_k.dt  <- tau.dt[amy.dt, on = .(PTID, VISIT)
                      ][triad.dt, on = .(PTID, VISIT)
                      ][!is.na(MOCA_score) & !is.na(AMY_1) & !is.na(TAU_1)]
rm(amy.dt, tau.dt)

## Model definitions
amy_ftrs    <- paste("AMY", 1:K_amy, sep = "_")
tau_ftrs    <- paste("TAU", 1:K_tau, sep = "_")

amy_covs    <- paste(amy_ftrs, collapse = " + ")
tau_covs    <- paste(tau_ftrs, collapse = " + ")

amy_reg     <- sprintf("%s ~ SEX_n + AGE_scan", amy_ftrs)
tau_reg     <- sprintf("%s ~ %s + SEX_n + AGE_scan", tau_ftrs, amy_covs)

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

hvr_kl.mod  <-
  str_glue("
# Latent variables
AMYLOID =~ {paste(amy_ftrs, collapse = ' + ')}
TAU =~ {paste(tau_ftrs, collapse = ' + ')}
{hvr.mod}
")

hvr_kd.mod  <-
  str_glue("
## Regressions ##
# Amyloid
{paste(amy_reg, collapse = '\n')}
# Tau
{paste(tau_reg, collapse = '\n')}
HVR_mean_inv ~ {amy_covs} + {tau_covs} + SEX_n + AGE_scan
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

moca_kl.mod <-
  str_glue("
# Latent variables
AMYLOID =~ {paste(amy_ftrs, collapse = ' + ')}
TAU =~ {paste(tau_ftrs, collapse = ' + ')}
{moca.mod}
")

moca_kd.mod <-
  str_glue("
## Regressions ##
# Amyloid
{paste(amy_reg, collapse = '\n')}
# Tau
{paste(tau_reg, collapse = '\n')}
HVR_mean_inv ~ {amy_covs} + {tau_covs} + SEX_n + AGE_scan
MOCA_score ~ {amy_covs} + {tau_covs} + HVR_mean_inv + SEX_n + AGE_scan + EDUC
")

## Fit models
# RMSEA = 0; CFI & TLI = 1 #
fname <- here("data/rds/mediation_selected_hvr.rds")
if (!file.exists(fname) | refit_mods) {
  mod_sel_hvr.fit <- sem(hvr.mod, data = triad_s.dt, #group = "DX",
                                estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(mod_sel_hvr.fit, fname)
} else {
  mod_sel_hvr.fit <- read_rds(fname)
}
rm(fname)

# RMSEA = 0; CFI & TLI = 1 #
fname <- here("data/rds/mediation_selected_moca.rds")
if (!file.exists(fname) | refit_mods) {
  mod_sel_moca.fit  <- sem(moca.mod, data = triad_s.dt, #group = "DX",
                                 estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(mod_sel_moca.fit, fname)
} else {
  mod_sel_moca.fit  <- read_rds(fname)
}
rm(fname)

# RMSEA > .18; low CFI/TLI
fname <- here("data/rds/mediation_latent_hvr.rds")
if (!file.exists(fname) | refit_mods) {
  mod_lat_hvr.fit <- sem(hvr_kl.mod, data = triad_k.dt, #group = "DX",
                                estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(mod_lat_hvr.fit, fname)
} else {
  mod_lat_hvr.fit <- read_rds(fname)
}
rm(fname)

# RMSEA > .16; low CFI/TLI
fname <- here("data/rds/mediation_latent_moca.rds")
if (!file.exists(fname) | refit_mods) {
  mod_lat_moca.fit <- sem(moca_kl.mod, data = triad_k.dt, #group = "DX",
                                 estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(mod_lat_moca.fit, fname)
} else {
  mod_lat_moca.fit <- read_rds(fname)
}
rm(fname)

# Unusable: RMSEA > .5; CFI < .3; TLI < 0
fname <- here("data/rds/mediation_clusters_hvr.rds")
if (!file.exists(fname) | refit_mods) {
  mod_clust_hvr.fit <- sem(hvr_kd.mod, data = triad_k.dt, #group = "DX",
                                estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(mod_clust_hvr.fit, fname)
} else {
  mod_clust_hvr.fit <- read_rds(fname)
}
rm(fname)

# Unusable: RMSEA ~ .5; CFI < .32; TLI < 0
fname <- here("data/rds/mediation_clusters_moca.rds")
if (!file.exists(fname) | refit_mods) {
  mod_clust_moca.fit <- sem(moca_kd.mod, data = triad_k.dt, #group = "DX",
                                 estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(mod_clust_moca.fit, fname)
} else {
  mod_clust_moca.fit <- read_rds(fname)
}
rm(fname)

### Simplest model with ALL ROIs
## HVR mediation
#fname <- here("data/rds/mediation_simplest_hvr.rds")
#if (!file.exists(fname) | refit_mods) {
  #model_simplest_hvr.fit <- sem(hvr.mod,
                                #data = triad_all.dt,
                                ##group = "DX",
                                #estimator = "ML")
                     ##se = "bootstrap", bootstrap = 10000)
  #write_rds(model_simplest_hvr.fit, fname)
  #semTable(model_simplest_hvr.fit, print.results = FALSE,
           #paramSets  = c("slopes", "constructed", "fits"),
           #columns    = c("estsestars", "z"),
           #fits       = c("npar", "chisq", "pvalue", "cfi", "tli",
                          #"aic", "bic2", "rmsea"),
           #type       = "html",
           #file       = here("data/derivatives/mediation_simplest_hvr.html"))
#} else {
  #model_simplest_hvr.fit <- read_rds(fname)
#}
#rm(fname)

## MoCA mediation
#fname <- here("data/rds/mediation_rois_simplest_moca.rds")
#if (!file.exists(fname) | refit_mods) {
  #model_simplest_moca.fit <- sem(moca.mod,
                                 #data = triad_all.dt,
                                 ##group = "DX",
                                 #estimator = "ML")
                     ##se = "bootstrap", bootstrap = 10000)
  #write_rds(model_simplest_moca.fit, fname)
  #semTable(model_simplest_moca.fit, print.results = FALSE,
           #paramSets  = c("slopes", "constructed", "fits"),
           #columns    = c("estsestars", "z"),
           #fits       = c("npar", "chisq", "pvalue", "cfi", "tli",
                          #"aic", "bic2", "rmsea"),
           #type       = "html",
           #file       = here("data/derivatives/mediation_simplest_moca.html"))
#} else {
  #model_simplest_moca.fit <- read_rds(fname)
#}
#rm(fname)

### Simple model with selected ROIs
## HVR mediation
#fname <- here("data/rds/mediation_simple_hvr.rds")
#if (!file.exists(fname) | refit_mods) {
  #model_simple_hvr.fit <- sem(hvr.mod,
                               #data = triad_sel.dt,
                               ##group = "DX",
                               #estimator = "ML")
                     ##se = "bootstrap", bootstrap = 10000)
  #write_rds(model_simple_hvr.fit, fname)
  #semTable(model_simple_hvr.fit, print.results = FALSE,
           #paramSets  = c("slopes", "constructed", "fits"),
           #columns    = c("estsestars", "z"),
           #fits       = c("npar", "chisq", "pvalue", "cfi", "tli",
                          #"aic", "bic2", "rmsea"),
           #type       = "html",
           #file       = here("data/derivatives/mediation_simple_hvr.html"))
#} else {
  #model_simple_hvr.fit <- read_rds(fname)
#}
#rm(fname)

## MoCA mediation
#fname <- here("data/rds/mediation_rois_simple_moca.rds")
#if (!file.exists(fname) | refit_mods) {
  #model_simple_moca.fit <- sem(moca.mod,
                               #data = triad_sel.dt,
                               ##group = "DX",
                               #estimator = "ML")
                     ##se = "bootstrap", bootstrap = 10000)
  #write_rds(model_simple_moca.fit, fname)
  #semTable(model_simple_moca.fit, print.results = FALSE,
           #paramSets  = c("slopes", "constructed", "fits"),
           #columns    = c("estsestars", "z"),
           #fits       = c("npar", "chisq", "pvalue", "cfi", "tli",
                          #"aic", "bic2", "rmsea"),
           #type       = "html",
           #file       = here("data/derivatives/mediation_simple_moca.html"))
#} else {
  #model_simple_moca.fit <- read_rds(fname)
#}
#rm(fname)

### Labels
#labels_cov  <- c(AGE_scan = "Age", SEX_n = "Sex", APOE_n = "APOE4")
#labels_hcv  <- c(HCv_l = "Left", HCv_r = "Right")
#labels_hvr  <- c(HVR = "HC-atrophy",
                 #HVR_lr = "1-HVR (Left)",
                 #HVR_rr = "1-HVR (Right)")
#labels_hvr2 <- c(HVR_mean_inv = "HC-atrophy")
#labels_moca <- c(MOCA_score = "MoCA")

### Models with latent variables ###
#amy_ftrs  <- paste("AMY", 1:5, sep = "_")
#tau_ftrs  <- paste("TAU", 1:5, sep = "_")

#latent.mod <-
  #str_glue("
           ## Latent variables
           #AMYLOID =~ {paste(amy_ftrs, collapse = ' + ')}
           #TAU =~ {paste(tau_ftrs, collapse = ' + ')}
           ## Regressions
           #AMYLOID ~ SEX_n + AGE_scan
           #TAU ~ a * AMYLOID + SEX_n + AGE_scan
           #HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan
           #MOCA_score ~ d * AMYLOID + e * TAU + f * HVR_mean_inv + SEX_n + AGE_scan + EDUC
           ## HVR mediation
           ## Direct effect
           #dAMY_HVR  := b
           ## Indirect effect
           #iTAU_HVR  := a * c
           ## Total effect
           #Total_HVR := dAMY_HVR + iTAU_HVR
           ## MOCA mediation
           ## Direct effect
           #dAMY_MOCA := d
           ## Indirect effects
           #iTAU_MOCA := a * e
           #iHVR_MOCA := Total_HVR * f
           ## Total effects
           #iTotal    := iTAU_MOCA + iHVR_MOCA
           #Total_MOCA := iTotal + dAMY_MOCA
           ## Proportion analysis
           #propAMY_HVR := dAMY_HVR / Total_HVR
           #propTAU_HVR := iTAU_HVR / Total_HVR
           #propAMY_MOCA := dAMY_MOCA / Total_MOCA
           #propTAU_MOCA := iTAU_MOCA / Total_MOCA
           #propHVR_MOCA := iHVR_MOCA / Total_MOCA
           #propIndirect_MOCA := iTotal / Total_MOCA
           #")

#fname <- here("data/rds/mediation_rois_latent.rds")
#if (!file.exists(fname) | refit_mods) {
  #model_rois_latent.fit <- sem(latent.mod,
                               #data = triad.dt[order(DX)],
                               ##group = "DX",
                               #estimator = "ML")
                     ##se = "bootstrap", bootstrap = 10000)
  #write_rds(model_rois_latent.fit, fname)
#} else {
  #model_rois_latent.fit <- read_rds(fname)
#}
#rm(fname)

#### Model with all mediators ###
## Text parsing
#amy_covars    <- paste(amy_ftrs, collapse = " + ")
#tau_covars    <- paste(tau_ftrs, collapse = " + ")

#amy_regress   <- sprintf("%s ~ SEX_n + AGE_scan", amy_ftrs)
#tau_regress   <- sprintf("%s ~ %s + SEX_n + AGE_scan",
                         #tau_ftrs, amy_covars)

## Model definition
#detailed.mod <-
  #str_glue("
### Regressions ##
## Amyloid
#{paste(amy_regress, collapse = '\n')}
## Tau
#{paste(tau_regress, collapse = '\n')}
#HVR_mean_inv ~ {amy_covars} + {tau_covars} + SEX_n + AGE_scan
#MOCA_score ~ {amy_covars} + {tau_covars} + HVR_mean_inv + SEX_n + AGE_scan + EDUC
#")

#fname <- here("data/rds/mediation_rois_detailed.rds")
#if (!file.exists(fname) | refit_mods) {
  #model_rois_detailed.fit <- sem(detailed.mod,
                                 #data = triad.dt[order(DX)],
                                 ##group = "DX",
                                 #estimator = "ML")
                     ##se = "bootstrap", bootstrap = 10000)
  #write_rds(model_rois_detailed.fit, fname)
#} else {
  #model_rois_detailed.fit <- read_rds(fname)
#}
#rm(fname)

## Estimates data.table
#estimates.dt  <- parameterEstimates(model_rois_detailed.fit,
                                    #standardized = TRUE) |>
           #as.data.table()

## Remove variances data
#estimates.dt  <- estimates.dt[op != "~~"]

## Adjust for multiple comparisons
#estimates.dt[, pvalue_adj := p.adjust(pvalue, method = "bonferroni")]

#estimates.dt[, lhs_roi_id := as.numeric(str_extract(lhs, "\\d{3}"))]
#estimates.dt[, rhs_roi_id := as.numeric(str_extract(rhs, "\\d{3}"))]
