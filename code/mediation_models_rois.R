#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
library(stringr)
library(lavaan)
library(lavaanPlot)
library(semTable)

## Refit models
#refit_models <- TRUE
refit_models  <- TRUE

## Print plots ### Needs to be done outside renv
print_plots   <- FALSE
#print_plots <- TRUE

## INPUT
# Load baseline data
fpath         <- here("data/rds/triad.rds")
if (file.exists(fpath)) {
  triad.dt    <- read_rds(fpath)
} else {
  here("code/demographics.R") |> source()
}
rm(fpath)

# Load Cerebra data
# Biomarkers
fpath         <- here("data/pet_biomarkers_cerebra.csv")
if (file.exists(fpath)) {
  cerebra.dt  <- fread(fpath)
} else {
  here("code/parse_csv_data.R") |> source()
}

# ROIs
fpaths        <- here("data/rds",
                      sprintf("cerebra_rois_%s_moca.rds", c("amy", "tau")))
if (any(!file.exists(fpaths))) {
  here("code/feature_selection.R") |> source()
} else {
  rois_amy    <- read_rds(fpaths[1])
  rois_tau    <- read_rds(fpaths[2])
}
rm(fpaths)

## Data cleaning
# Convert Sex to dummy variable
triad.dt[, SEX_n := as.numeric(SEX) - 1]

# Remove youth and dementias
triad.dt      <- triad.dt[!DX_clean %in% c("Young", "Other", "AD")]

# 1 - HVR (average for both sides)
#triad.dt[, `:=`(HVR_lr = 1 - HVR_l, HVR_rr = 1 - HVR_r)]
triad.dt[, `:=`(HVR_mean_inv = 1 - (HVR_l + HVR_r) / 2)]

triad.dt      <- triad.dt[, .(PTID, VISIT, DX, SEX_n, AGE_scan, EDUC,#APOE_n,
                              HVR_mean_inv, MOCA_score)]

# Merge triad and cerebra data
# Cerebra dictionary
dict_roi      <- unique(cerebra.dt[, .(LABEL_id, LABEL_name, SIDE)])

# Amyloid weighted average of all regions
amy.dt        <- cerebra.dt[!is.na(AMYLOID),
                            .(AMYLOID = weighted.mean(AMYLOID, VOL)),
                            .(PTID, VISIT)]
amy.dt        <- amy.dt[triad.dt, on = .(PTID, VISIT)]

# Tau weighted average of all regions
tau.dt        <- cerebra.dt[!is.na(TAU),
                            .(TAU = weighted.mean(TAU, VOL)),
                            .(PTID, VISIT)]

triad_all.dt  <- tau.dt[amy.dt, on = .(PTID, VISIT)
                        ][!is.na(MOCA_score) & !is.na(AMYLOID) & !is.na(TAU)]
rm(amy.dt, tau.dt)

# Amyloid weighted average of selected features
#amy.dt      <- rois_amy.dt[!is.na(cluster),
                           #.(LABEL_id = as.numeric(str_extract(id, "\\d{3}")),
                             #CLUSTER  = sprintf("AMY_%i", cluster))
                           #][cerebra.dt, on = "LABEL_id"
                           #][!is.na(AMYLOID) & !is.na(CLUSTER),
                           #weighted.mean(AMYLOID, VOL),
                           #.(PTID, VISIT, CLUSTER)] |>
  #dcast(... ~ CLUSTER, value.var = "V1")
amy.dt        <-
  cerebra.dt[rois_amy.dt[decision == "Confirmed",
                         .(LABEL_id = as.numeric(str_extract(id, "\\d{3}")))],
             on = "LABEL_id"
             ][!is.na(AMYLOID),
               .(AMYLOID = weighted.mean(AMYLOID, VOL)),
               .(PTID, VISIT)]

amy.dt        <- amy.dt[triad.dt, on = .(PTID, VISIT)]

# Tau weighted average of selected features
#tau.dt      <- rois_tau.dt[!is.na(cluster),
                           #.(LABEL_id = as.numeric(str_extract(id, "\\d{3}")),
                             #CLUSTER  = sprintf("TAU_%i", cluster))
                           #][cerebra.dt, on = "LABEL_id"
                           #][!is.na(TAU) & !is.na(CLUSTER),
                           #weighted.mean(TAU, VOL),
                           #.(PTID, VISIT, CLUSTER)] |>
  #dcast(... ~ CLUSTER, value.var = "V1")
tau.dt        <-
  cerebra.dt[rois_tau.dt[decision == "Confirmed",
                         .(LABEL_id = as.numeric(str_extract(id, "\\d{3}")))],
             on = "LABEL_id"
             ][!is.na(TAU),
               .(TAU = weighted.mean(TAU, VOL)),
               .(PTID, VISIT)]

triad_sel.dt  <- tau.dt[amy.dt, on = .(PTID, VISIT)
                        ][!is.na(MOCA_score) & !is.na(AMYLOID) & !is.na(TAU)]
#triad.dt  <- triad.dt[!is.na(AMY_1) & !is.na(TAU_1) & !is.na(MOCA_score)]
rm(amy.dt, tau.dt)

## Simple model with selected ROIs
## Model definitions
hvr.mod       <-
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

moca.mod <-
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

## Simplest model with ALL ROIs
# HVR mediation
fname <- here("data/rds/mediation_simplest_hvr.rds")
if (!file.exists(fname) | refit_models) {
  model_simplest_hvr.fit <- sem(hvr.mod,
                                data = triad_all.dt,
                                #group = "DX",
                                estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(model_simplest_hvr.fit, fname)
  semTable(model_simplest_hvr.fit, print.results = FALSE,
           paramSets  = c("slopes", "constructed", "fits"),
           columns    = c("estsestars", "z"),
           fits       = c("npar", "chisq", "pvalue", "cfi", "tli",
                          "aic", "bic2", "rmsea"),
           type       = "html",
           file       = here("data/derivatives/mediation_simplest_hvr.html"))
} else {
  model_simplest_hvr.fit <- read_rds(fname)
}
rm(fname)

# MoCA mediation
fname <- here("data/rds/mediation_rois_simplest_moca.rds")
if (!file.exists(fname) | refit_models) {
  model_simplest_moca.fit <- sem(moca.mod,
                                 data = triad_all.dt,
                                 #group = "DX",
                                 estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(model_simplest_moca.fit, fname)
  semTable(model_simplest_moca.fit, print.results = FALSE,
           paramSets  = c("slopes", "constructed", "fits"),
           columns    = c("estsestars", "z"),
           fits       = c("npar", "chisq", "pvalue", "cfi", "tli",
                          "aic", "bic2", "rmsea"),
           type       = "html",
           file       = here("data/derivatives/mediation_simplest_moca.html"))
} else {
  model_simplest_moca.fit <- read_rds(fname)
}
rm(fname)

## Simple model with selected ROIs
# HVR mediation
fname <- here("data/rds/mediation_simple_hvr.rds")
if (!file.exists(fname) | refit_models) {
  model_simple_hvr.fit <- sem(hvr.mod,
                               data = triad_sel.dt,
                               #group = "DX",
                               estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(model_simple_hvr.fit, fname)
  semTable(model_simple_hvr.fit, print.results = FALSE,
           paramSets  = c("slopes", "constructed", "fits"),
           columns    = c("estsestars", "z"),
           fits       = c("npar", "chisq", "pvalue", "cfi", "tli",
                          "aic", "bic2", "rmsea"),
           type       = "html",
           file       = here("data/derivatives/mediation_simple_hvr.html"))
} else {
  model_simple_hvr.fit <- read_rds(fname)
}
rm(fname)

# MoCA mediation
fname <- here("data/rds/mediation_rois_simple_moca.rds")
if (!file.exists(fname) | refit_models) {
  model_simple_moca.fit <- sem(moca.mod,
                               data = triad_sel.dt,
                               #group = "DX",
                               estimator = "ML")
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(model_simple_moca.fit, fname)
  semTable(model_simple_moca.fit, print.results = FALSE,
           paramSets  = c("slopes", "constructed", "fits"),
           columns    = c("estsestars", "z"),
           fits       = c("npar", "chisq", "pvalue", "cfi", "tli",
                          "aic", "bic2", "rmsea"),
           type       = "html",
           file       = here("data/derivatives/mediation_simple_moca.html"))
} else {
  model_simple_moca.fit <- read_rds(fname)
}
rm(fname)
