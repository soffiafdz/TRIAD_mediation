#!/usr/bin/env Rscript

library(here)
library(data.table)
library(stringr)
library(lubridate)
library(progress) ## needed?


### CONSTANTS
USE_IMPUTED <- T

### INPUTS
fpaths      <- here("data", c(
  "derivatives/adni_icc_scale.csv",
  "rds/adni_dxs_imputed.rds",
  "rds/adni_vols-hcvc.rds"
))

## ICC/Scale
if (!file.exists(fpaths[1])) {
  sprintf("File: %s is required but could not be found.", fpaths) |> stop()
}

scale.dt    <- fpaths[1] |> fread() |> setnames("VISIT", "EXAMDATE")
scale.dt[, EXAMDATE := ymd(EXAMDATE)]
setkey(scale.dt, PTID, EXAMDATE)

## Imputed Dxs from ADNIMERGE
if (file.exists(fpaths[2])) {
  dx.dt     <- readRDS(fpaths[2])
} else {
  here("code/impute_dx.R") |> source()
}

## HcVc volumes (QCed)
if (file.exists(fpaths[3])) {
  hcvc.dt   <- readRDS(fpaths[3])
} else {
  here("code/qc-filter_segms_adni.R") |> source()
}

rm(fpaths)

### Data PROCESSING
## Bring back volumes to native CC
## There is a subject with SCALE 0; remove
vols.dt     <- scale.dt[
  hcvc.dt
][
  SCALEFACTOR != 0 & ICC != 0,
  .(
    PTID,
    EXAMDATE,
    ICC   = ICC  / 1000,
    HC_l  = LHC  / (SCALEFACTOR * 1000),
    HC_r  = RHC  / (SCALEFACTOR * 1000),
    VC_l  = LCSF / (SCALEFACTOR * 1000),
    VC_r  = RCSF / (SCALEFACTOR * 1000)
  )
]
#rm(hcvc.dt, scale.dt)

## Use imputed Dx?
if (USE_IMPUTED) {
  dx.dt[, DX := str_remove(DX, "\\?")]
} else {
  dx.dt <- dx.dt[!DX %like% "\\?"]
}

vols.dt     <- vols.dt |>
  melt(measure = patterns("HC|VC"), value = "CC") |>
  {
    function(DT)
    DT[
      , c("ROI", "SIDE") := tstrsplit(variable, split = "_")
    ][
      , SIDE := toupper(SIDE)
    ][
      , variable := NULL
    ]
  }() |>
  setkey(PTID, EXAMDATE) |>
  merge(dx.dt) |>
  setcolorder("CC", after = "SIDE") |>
  setcolorder(c("VISCODE", "DX"), before = "ICC")

### Head-size adjustment methods
hc_hvr.dt <- vols.dt[
  , B := summary(lm(CC ~ ICC))[[4]][2], .(ROI, SIDE)
][
  , let(
    NON_cc = CC,
    PRP_cc = CC / ICC,
    RES_cc = CC - B * (ICC - vols.dt["CN", on = "DX", mean(ICC)])
  )
][
  , c("CC", "B") := NULL
] |>
  melt(measure = patterns("_cc$"), variable = "ADJ") |>
  dcast(... ~ ROI, value = "value") |>
  # HVR
  {
    function(DT)
    DT[
      , ADJ := str_remove(ADJ, "_cc$")
    ][
      !"NON", on = "ADJ", HVR := (HC / (HC + VC))
    ]
  }() |>
  setkey(PTID, EXAMDATE)

### OUTPUT
outrds      <- here("data/rds/adni_hc-hvr_adj.rds")
saveRDS(hc_hvr.dt, outrds)
rm(outrds)
