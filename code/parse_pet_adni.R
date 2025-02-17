#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lubridate)


### INPUTS
petdir      <- here("data/source_adni_pet_data/csv_files")
fpaths      <- c(here("data/rds/adni_hc-hvr.rds"),
                 list.files(petdir, pattern = "BERKELEY", full = TRUE))
rm(petdir)

for (fpath in fpaths) {
  if (!file.exists(fpath)) {
    sprintf("File: %s is required but could not be found.", fpath) |> stop()
  }
  rm(fpath)
}

hcvr.dt     <- readRDS(fpaths[1])
ucb_amy.dt  <- fread(fpaths[2], key = "PTID", select = c(5:6,9:10,12,14))
ucb_tau.dt  <- fread(fpaths[3], key = "PTID", select = c(5:6,9,12))
#ucb_tau2.dt <- fread(fpaths[4]) ## Partial volume correction
#ucb_fdg.dt <- fread(fpaths[5])
rm(fpaths)

## Remove faulty scans
ucb_amy.dt  <- ucb_amy.dt[qc_flag > 0, -"qc_flag"]
ucb_tau.dt  <- ucb_tau.dt[qc_flag > 0, -"qc_flag"]

## Merge on closest date
dts_amy.dt  <- ucb_amy.dt[, .(DATE_AMY = ymd(SCANDATE)), keyby = "PTID"]
dts_tau.dt  <- ucb_tau.dt[, .(DATE_TAU = ymd(SCANDATE)), keyby = "PTID"]

## tau < amy; merge ON amy
dts1.dt     <- dts_tau.dt[dts_amy.dt, allow.cartesian = TRUE] |> na.omit()
dts1.dt[, DIFF_tau_amy := DATE_AMY - DATE_TAU]
dts1.dt     <- dts1.dt[, .SD[which.min(abs(DIFF_tau_amy))], .(PTID, DATE_TAU)
                       ][, .SD[which.min(abs(DIFF_tau_amy))],
                       .(PTID, DATE_AMY)]
ucb_tau.dt  <- dts1.dt[ucb_tau.dt, on = .(PTID, DATE_TAU = SCANDATE)]

## amy < mri; merge on mri
dts_mri.dt  <- hcvr.dt[, .(DATE_MRI = ymd(EXAMDATE)), keyby = "PTID"]
dts2.dt     <- dts_amy.dt[dts_mri.dt, allow.cartesian = TRUE] |> na.omit()
dts2.dt[, DIFF_amy_mri := DATE_AMY - DATE_MRI]
dts2.dt     <- dts2.dt[, .SD[which.min(abs(DIFF_amy_mri))], .(PTID, DATE_AMY)
                       ][, .SD[which.min(abs(DIFF_amy_mri))],
                       .(PTID, DATE_MRI)]
ucb_amy.dt  <- dts2.dt[ucb_amy.dt, on = .(PTID, DATE_AMY = SCANDATE)]

## tau < mri; merge on mri
dts3.dt     <- dts_tau.dt[dts_mri.dt, allow.cartesian = TRUE] |> na.omit()
dts3.dt[, DIFF_tau_mri := DATE_TAU - DATE_MRI]
dts3.dt     <- dts3.dt[, .SD[which.min(abs(DIFF_tau_mri))], .(PTID, DATE_TAU)
                       ][, .SD[which.min(abs(DIFF_tau_mri))],
                       .(PTID, DATE_MRI)]
ucb_tau.dt  <- dts3.dt[ucb_tau.dt, on = .(PTID, DATE_TAU)]
#rm(dts1.dt, dts2.dt, dts3.dt, dts_amy.dt, dts_tau.dt, dts_mri.dt, hcvr.dt)

## Output
outrds      <- here('data/rds', sprintf("ucb_pet-%s.rds", c("amy", "tau")))
for (i in 1:2) readr::write_rds(list(ucb_amy.dt, ucb_tau.dt)[[i]], outrds[i])
#rm(i, outrds)
