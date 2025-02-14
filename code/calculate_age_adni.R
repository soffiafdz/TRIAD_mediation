#!/usr/bin/env Rscript

library(here)
library(data.table)
library(stringr)
library(lubridate)
library(ADNIMERGE)


### INPUT
fpaths      <- here("data", c("/derivatives/vols_hcvc_adni.csv",
                              "MRIMETA_27Aug2024.csv",
                              "MRI3META_27Aug2024.csv"))
for (fpath in fpaths) {
  if (!file.exists(fpath)) {
    sprintf("File: %s is required but could not be found.", fpath) |> stop()
  }
  rm(fpath)
}

hcvc.dt     <- fread(fpaths[1])
mrimeta.dt  <- fread(fpaths[2])
mrimeta3.dt <- fread(fpaths[3])
rm(fpaths)

## Obtain AGE from ADNIMERGE
data(adnimerge)
setDT(adnimerge)
age_bl.dt   <- adnimerge[!is.na(AGE)
                         ][VISCODE == "bl",
                         .(PTID, EXAMDATE, AGE)] |>
  setnames(c("EXAMDATE", "AGE"), c("EXAMDATE_bl", "AGE_bl")) |>
  setkey(PTID)
rm(adnimerge)

## Concatenate mrimeta files
mrimeta.dt  <- rbind(mrimeta.dt[, .(PTID, EXAMDATE, VISCODE2)],
                     mrimeta3.dt[, .(PTID, EXAMDATE, VISCODE2)]) |>
  unique() |> setkey(PTID, EXAMDATE)
rm(mrimeta3.dt)

## Extract PTID and EXAMDATE from the filename of hcvc.dt
mrisubs.dt  <-
  hcvc.dt[, .(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
              EXAMDATE = ymd(str_extract(ID, "\\d{4}(-\\d{2}){2}|\\d{8}")))] |>
  setkey(PTID, EXAMDATE)
rm(hcvc.dt)

### MERGE data.tables
age.dt      <- age_bl.dt[mrisubs.dt] |> setkey(PTID, EXAMDATE)
age.dt      <- mrimeta.dt[age.dt]
rm(age_bl.dt, mrimeta.dt, mrisubs.dt)

## Remove duplicate rows from merging
# 070_S_5040: 2013-07-25 is both m03 & m06; yet there is already a m03
# 130_S_0505: 2006-10-18 is both sc & bl; yet there is another, older, sc
age.dt      <- age.dt[!((PTID == "070_S_5040"
                         & VISCODE2 == "m03"
                         & EXAMDATE == "2013-07-25") |
                        (PTID == "130_S_0505"
                         & VISCODE2 == "sc"
                         & EXAMDATE == "2006-10-18"))]

### CALCULATE AGE using dates (from MRI sessions)
age.dt[, `:=`(EXAMDATE = ymd(EXAMDATE), EXAMDATE_bl = ymd(EXAMDATE_bl))]
age.dt[, DIFF_years := as.numeric(EXAMDATE - EXAMDATE_bl) / 365.4]
age.dt[, AGE := AGE_bl + DIFF_years]
age.dt |> setnames("VISCODE2", "VISCODE") |> setkey(PTID, EXAMDATE)


### FILTER MRI sub

## MERGE data.tables
age.dt      <- age.dt[, .(PTID, VISCODE, EXAMDATE, AGE)]

### OUTPUT
outrds      <- here("data/rds/adni_age_calculated.rds")
age.dt |> readr::write_rds(outrds)
rm(outrds)
