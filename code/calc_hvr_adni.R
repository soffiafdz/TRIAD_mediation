#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(lubridate)
library(stringr)


### INPUTS
fpaths      <- here("data", c("derivatives/icc_scale_adni_old.csv",
                              "rds/adni_dxs_imputed.rds",
                              "rds/adni_vols-hcvc.rds"))

## ICC/Scale
if (!file.exists(fpaths[1])) {
  sprintf("File: %s is required but could not be found.", fpaths) |> stop()
}

scale.dt    <- fpaths[1] |> fread() |> setnames("VISIT", "EXAMDATE")
scale.dt[, EXAMDATE := ymd(EXAMDATE)]
setkey(scale.dt, PTID, EXAMDATE)

## Imputed Dxs from ADNIMERGE
if (file.exists(fpaths[2])) {
  dx.dt     <- read_rds(fpaths[2])
} else {
  here("code/impute_dx.R") |> source()
}

## HcVc volumes (QCed)
if (file.exists(fpaths[3])) {
  hcvc.dt   <- read_rds(fpaths[3])
} else {
  here("code/qc-filter_segms_adni.R") |> source()
}

rm(fpaths)

### Data PROCESSING
## Bring back volumes to native CC
## There is a subject with SCALE 0; remove
vols.dt     <- scale.dt[hcvc.dt
                        ][SCALEFACTOR != 0 & ICC != 0,
                        .(PTID, EXAMDATE,
                          ICC    = ICC  / 1000,
                          HC_l   = LHC  / (SCALEFACTOR * 1000),
                          HC_r   = RHC  / (SCALEFACTOR * 1000),
                          CSF_l  = LCSF / (SCALEFACTOR * 1000),
                          CSF_r  = RCSF / (SCALEFACTOR * 1000))]
rm(hcvc.dt, scale.dt)

vols.dt     <- vols.dt |>
  melt(id.vars = c("PTID", "EXAMDATE", "ICC"),
       variable.name = "ROI", value.name = "CC") |>
  setkey(PTID, EXAMDATE)

## Adjust for ICV
## Linear model of VAL ~ ROI & average ICC on healthy people
## Confirmer or Imputed??
cn.dt       <- dx.dt[vols.dt
                     ][DX %like% "CN",
                     .SD[which.min(EXAMDATE)],
                     .(PTID, ROI)]
icc_cn      <- cn.dt[!duplicated(PTID), mean(ICC)]
b.dt        <- cn.dt[, .(b = summary(lm(CC ~ ICC))$coefficients[2]), ROI]

## Adjust by head size
vols.dt[, c("ROI", "SIDE") := tstrsplit(toupper(ROI), "_")]
vols.dt[, CC_adj := CC - b.dt[ROI == ROI, b] * (ICC - icc_cn)]
hc_hvr.dt   <- vols.dt[, -"ICC"] |>
  dcast(... ~ ROI, value.var = c("CC", "CC_adj"))
hc_hvr.dt[, HVR := CC_HC / (CC_HC + CC_CSF)]
hc_hvr.dt[, c("CC_CSF", "CC_HC", "CC_adj_CSF") := NULL]
setnames(hc_hvr.dt, "CC_adj_HC", "HCvol_adj")
rm(b.dt, icc_cn, vols.dt)


### OUTPUT
outrds      <- here("data/rds/adni_hc-hvr.rds")
write_rds(hc_hvr.dt, outrds)
rm(outrds)
