#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lubridate)
library(readr)

### INPUT
fpaths      <- here("data", c("derivatives/vols_hcvc_adni.csv",
                              "qc_adni_hcvc_outliers.csv",
                              "qc_adni_hcvc_outliers_dx.csv"))

if (any(!file.exists(fpaths))) {
  sprintf("File: %s is required but could not be found.", fpath) |> stop()
}

## Volumes
hcvc.dt     <- fread(fpaths[1])
hcvc.dt[, PTID := str_extract(ID, "\\d{3}_S_\\d{4}")]
hcvc.dt[, EXAMDATE := ymd(str_extract(ID, "\\d{4}(-\\d{2}){2}|\\d{8}"))]
hcvc.dt[, ID := NULL]
setkey(hcvc.dt, PTID, EXAMDATE)

## QC list
qc.dt       <- fpaths[-1] |>
  lapply(fread, header = FALSE) |>
  rbindlist() |>
  setnames(c("ID", "QC"))

qc.dt[, PTID := str_extract(ID, "\\d{3}_S_\\d{4}")]
qc.dt[, EXAMDATE := ymd(str_extract(ID, "\\d{4}(-\\d{2}){2}|\\d{8}"))]
qc.dt[, ID := NULL]
setkey(qc.dt, PTID, EXAMDATE)


### QC filtering
## Warning & Fail
hcvc.dt     <- hcvc.dt[!qc.dt[QC %in% c("Warning", "Fail")]] |>
  setcolorder(7:8)


### OUTPUT
outrds      <- here("data/rds/adni_vols-hcvc.rds")
readr::write_rds(hcvc.dt, outrds)
rm(fpaths, qc.dt, outrds)
