#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
library(lavaan)
library(parameters)

## Imputed Dxs from ADNIMERGE
REFIT_FA    <- T

## PET
fpaths      <- here("data/rds", sprintf("ucb_pet-%s.rds", c("amy", "tau")))
if (all(file.exists(fpaths))) {
  ucb_amy.dt <- readRDS(fpaths[1])
  ucb_tau.dt <- readRDS(fpaths[2])
} else {
  here("code/parse_pet_adni.R") |> source()
}

## WMH
fpath       <- here("data/adni_wmh-vol.csv")
if (!file.exists(fpath)) {
  sprintf("File: %s is required but could not be found.", fpath) |> stop()
}
wmh.dt      <- fread(fpath)
rm(fpath)

## HCvol & HVR
fpath       <- here("data/rds/adni_hc-hvr.rds")
if (file.exists(fpath)) {
  hc_hvr.dt <- read_rds(fpath)
} else {
  here("code/calc_hvr_adni.R") |> source()
}

## Age (calculated)
fpath       <- here("data/rds/adni_age_calculated.rds")
if (file.exists(fpath)) {
  age.dt    <- read_rds(fpath)
} else {
  here("code/calculate_age_adni.R") |> source()
}

# Merge data
mri.dt      <- hc_hvr.dt[, .(HC = mean(HCvol_adj), HVR = mean(HVR)),
                         keyby = .(PTID, EXAMDATE)
                         ][wmh.dt[, .(WMHlog = log(WMHvol + 1)),
                                  keyby = .(PTID, EXAMDATE)]] |>
na.omit()

mri.dt      <- age.dt[, -"VISCODE"][mri.dt, on = .(PTID, EXAMDATE)] |>
na.omit()

pet.dt      <- ucb_amy.dt[ucb_tau.dt, on = .(PTID, DATE_AMY, DATE_MRI)] |>
na.omit()

pet.dt      <- age.dt[, -"VISCODE"][pet.dt, on = .(PTID, EXAMDATE = DATE_MRI)]


### Remove rows with higher time difference of 1y between scans
biomker.dt  <- mri.dt[pet.dt, on = .(PTID, EXAMDATE)
                      ][DIFF_tau_mri < 365 & DIFF_tau_amy < 365] |>
na.omit()

### CFA
## MRI
# Good fit (CFI, TLI, RMSE); but variance is negative.

cfa_mri.m   <- '
Neurodeg. =~ AGE + HVR + WMHlog
'

fpath <- here("data/rds/adni_mcfa-mri1.rds")
mri.dt[, AGE_scl := scale(AGE, scale = T)]
if (!file.exists(fpath) | REFIT_FA) {
  cfa1.f    <- cfa(mri.dt, model = cfa_mri.m, std.lv = T,
                   cluster = "PTID") |>
  write_rds(fpath)
} else {
  cfa1.f    <- read_rds(fpath)
}

## BIOMARKERS
# Good fit; CFI & TLI == 1.; RMSEA == 0
# High variances in Amyloid, but it's normal due to scale.
# Fitting the variance to 0 fixes the issue.
# WMH estimate is not significant
cfa_bio.m   <- '
Pathology =~ WMHlog + AMY_centiloids + TAU_metaROI
TAU_metaROI ~~ 0*TAU_metaROI
'

setnames(biomker.dt, c("META_TEMPORAL_SUVR", "CENTILOIDS"),
         c("TAU_metaROI", "AMY_centiloids"))
fpath <- here("data/rds/adni_mcfa-biomarkers1.rds")
if (!file.exists(fpath) | REFIT_FA) {
  cfa2.f    <- cfa(biomker.dt, model = cfa_bio.m, std.lv = T,
                   cluster = "PTID") |>
  write_rds(fpath)
} else {
  cfa2.f    <- read_rds(fpath)
}

outmd <- here('data/derivatives/adni_mcfa-biomarkers.md')
sink(outmd)
parameters(cfa2.f) |> print_md()
sink()
rm(outmd)

## Not a good fit: RMSEA ~.1
cfa_bio2.m   <- '
Pathology =~ AMY_centiloids + TAU_metaROI + HVR
TAU_metaROI ~~ 0*TAU_metaROI
'

fpath <- here("data/rds/adni_mcfa-pet_hvr.rds")
if (!file.exists(fpath) | REFIT_FA) {
  cfa2.2.f    <- cfa(biomker.dt, model = cfa_bio2.m, std.lv = T,
                   cluster = "PTID") |>
  write_rds(fpath)
} else {
  cfa2.2.f    <- read_rds(fpath)
}

outmd <- here('data/derivatives/adni_mcfa-pet_hvr.md')
sink(outmd)
parameters(cfa2.2.f) |> print_md()
sink()
rm(outmd)

cfa_bio3.m   <- '
Pathology =~ WMHlog + AMY_centiloids + TAU_metaROI + HVR
TAU_metaROI ~~ 0*TAU_metaROI
'

fpath <- here("data/rds/adni_mcfa-biomarkers_all.rds")
if (!file.exists(fpath) | REFIT_FA) {
  cfa2.3.f    <- cfa(biomker.dt, model = cfa_bio3.m, std.lv = T,
                   cluster = "PTID") |>
  write_rds(fpath)
} else {
  cfa2.3.f    <- read_rds(fpath)
}

outmd <- here('data/derivatives/adni_mcfa-biomarkers_all.md')
sink(outmd)
parameters(cfa2.3.f) |> print_md()
sink()
rm(outmd)
### Two factors:
### Good fit
### MCFA remains good; WMH estimate is not significant in MRI
cfa2_hvr.m  <- '
PET =~ AMY_centiloids + TAU_metaROI
MRI =~ HVR + WMHlog
'
fpath <- here("data/rds/adni_mcfa-biomarkers2_hvr.rds")
if (!file.exists(fpath) | REFIT_FA) {
  cfa3.1.f   <- cfa(biomker.dt, model = cfa2_hvr.m, std.lv = T,
                    cluster = "PTID") |>
  write_rds(fpath)
} else {
  cfa3.1.f     <- read_rds(fpath)
}

outmd <- here('data/derivatives/adni_mcfa-pet_mri-hvr.md')
sink(outmd)
parameters(cfa3.1.f) |> print_md()
sink()
rm(outmd)

## Using HC has a worse fit (Higher AIC * BIC)
fpath <- here("data/rds/adni_mcfa-biomarkers2_hc.rds")
cfa2_hc.m <- sub("HVR", "HC", cfa2_hvr.m)
if (!file.exists(fpath) | REFIT_FA) {
  cfa3.2.f   <- cfa(biomker.dt, model = cfa2_hc.m,
                    cluster = "PTID") |>
  write_rds(fpath)
} else {
  cfa3.2.f     <- read_rds(fpath)
}

outmd <- here('data/derivatives/adni_mcfa-pet_mri-hc.md')
sink(outmd)
parameters(cfa3.2.f) |> print_md()
sink()
rm(outmd)

## Second order factor
### Model not identified
#cfa3.m      <- '
#PET =~ AMY_centiloids + TAU_metaROI
#MRI =~ HVR + WMHlog
#Biomarkers =~ 1*PET + 1*MRI
#Biomarkers ~~ Biomarkers
#'

#fpath <- here("data/rds/adni_mcfa-biomarkers3_hc.rds")
#if (!file.exists(fpath) | REFIT_FA) {
  #cfa4.f   <- cfa(biomker.dt, model = cfa3.m, cluster = "PTID") |>
  #write_rds(fpath)
#} else {
  #cfa4.f     <- read_rds(fpath)
#}

### Save latent factors
path.dt  <- biomker.dt[, 1:2]
path.dt[, PATHOLOGY_lat := lavPredict(cfa2.f)]
path.dt[, PET_lat := lavPredict(cfa3.1.f)[,1]]
path.dt[, MRI_lat := lavPredict(cfa3.1.f)[,2]]

here('data/rds/adni_biomarkers-latent.rds') |>
readr::write_rds(x = path.dt)

