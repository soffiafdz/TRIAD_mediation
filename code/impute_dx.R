#!/usr/bin/env Rscript

library(here)
library(data.table)
library(stringr)
library(lubridate)
library(ADNIMERGE)
library(progress)

### INPUT
fpaths      <- here("data", c(
  "/derivatives/vols_hcvc_adni.csv",
  "MRIMETA_27Aug2024.csv",
  "MRI3META_27Aug2024.csv"
))

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

### CLEAN
## Extract PTID and EXAMDATE from the filename
mrisubs.dt  <-
  hcvc.dt[, .(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
              EXAMDATE = ymd(str_extract(ID, "\\d{4}(-\\d{2}){2}|\\d{8}")))] |>
  setkey(PTID, EXAMDATE)
rm(hcvc.dt)

## Obtain DX from ADNIMERGE
data(adnimerge)
setDT(adnimerge)
dx.dt       <- adnimerge[!is.na(DX), .(PTID, VISCODE, DX)] |>
  setkey(PTID, VISCODE)

## EXAMDATE from MRIMETA is needed for joining with ADNIMERGE
## Concatenate mrimeta files
mrimeta.dt  <- rbind(mrimeta.dt[, .(PTID, EXAMDATE, VISCODE2)],
                     mrimeta3.dt[, .(PTID, EXAMDATE, VISCODE2)]) |>
  setkey(PTID, VISCODE2)

dx.dt       <- dx.dt[mrimeta.dt] |> setkey(PTID, EXAMDATE)
dx.dt       <- dx.dt[mrisubs.dt] |> unique()
dx.dt[, DX := as.character(DX)]
rm(adnimerge, mrimeta.dt, mrimeta3.dt, mrisubs.dt)

## Remove duplicate rows from merging
# 070_S_5040: 2013-07-25 is both m03 & m06; yet there is already a m03
# 130_S_0505: 2006-10-18 is both sc & bl; yet there is another, older, sc
dx.dt       <- dx.dt[
  !(
    (PTID == "070_S_5040" & VISCODE == "m03" & EXAMDATE == "2013-07-25") |
    (PTID == "130_S_0505" & VISCODE == "sc" & EXAMDATE == "2006-10-18")
  )
]

### IMPUTATION ALGORITHM
## Progress bar
pb <- progress_bar$new(
  format = "Imputation | :what [:bar] :current/:total",
  total = dx.dt[is.na(DX), .N],
  clear = FALSE,
  width = 75
)

# Iterate through subjects with missing Dx
for (subj in dx.dt[is.na(DX), PTID]) {
  # Tick progress bar
  pb$tick(tokens = list(what = sprintf(
    "PTID: %s — missing %i Dxs",
    subj,
    dx.dt[PTID == subj][is.na(DX), .N]
  )))

  # Calculate the number of available Dx for i subject
  dx_n.dt   <- dx.dt[PTID == subj][!is.na(DX), .N, DX]

  # If there is not any usable Dx, skip subject
  if (dx_n.dt[, .N == 0]) next

  # If there is only one usable Dx, use that one
  if (dx_n.dt[, .N == 1]) {
    imp_dx  <- dx_n.dt[, DX]
    dx.dt[PTID == subj & is.na(DX), DX := paste0(imp_dx, "?")]
    next
  }

  # If there is more than one, compare with previous&following Dxs
  # Extract subsample of subject with a column for chronological order
  dt <- dx.dt[PTID == subj]
  dt[order(EXAMDATE), I := rowid(PTID)]

  # Reiterate algorithm until there are no more viable imputations
  while (dt[is.na(DX), .N != 0]) {
    # Iterate through visits with missing Dx
    for (i in dt[is.na(DX)][order(I), I]) {
      # Log visit and date
      visc_dx <- dt[I == i, VISCODE]
      date_dx <- dt[I == i, EXAMDATE]

      # Extract previous&following Dxs
      prev_dx <- if (i == 1) NA else dt[I == i - 1, DX]
      next_dx <- if (i == dt[, .N]) NA else dt[I == i + 1, DX]
      dxs     <- c(prev_dx, next_dx)

      # Check if both surrounding Dxs are missing
      if (all(is.na(dxs))) {
        # If viable, skip visit until other visits are imputed;
        # Else, skip subject altogether
        if (dt[is.na(DX), .N == 1]) break else next
      } else if (is.na(prev_dx)) {
        # If previous Dx is missing, impute with following Dx
        imp_dx <- fifelse(next_dx %like% "\\?", next_dx, paste0(next_dx, "?"))
      } else if (is.na(next_dx)) {
        # If following Dx is missing, impute with previous Dx
        imp_dx <- fifelse(prev_dx %like% "\\?", prev_dx, paste0(prev_dx, "?"))
      } else if (prev_dx == next_dx) {
        # If the Dxs are equal, impute with previous (no difference)
        imp_dx <- fifelse(prev_dx %like% "\\?", prev_dx, paste0(prev_dx, "?"))
      } else if (sum(dxs %like% "\\?") == 1) {
        # If only one of the Dx is itself imputed, impute with the other one
        imp_dx <- fifelse(prev_dx %like% "\\?", next_dx, prev_dx)
      } else {
        # If unequal Dxs; whether imputed or not, use the Dx closer in date
        imp_dx <- dt[I %in% c(i-1, i+1), .(DX, DIFF = EXAMDATE - date_dx)
                     ][order(abs(DIFF)), DX][1]
        if (!imp_dx %like% "\\?") imp_dx <- paste0(imp_dx, "?")
      }

      dt[I == i, DX := imp_dx]
      dx.dt[PTID == subj & VISCODE == visc_dx, DX := imp_dx]
      rm(date_dx, dxs, i, next_dx, prev_dx, visc_dx)
    }
  }
}
dx.dt   <- dx.dt[!is.na(DX)]
rm(dt, dx_n.dt, imp_dx, subj, pb)

### OUTPUT
## RDS
outrds  <- here("data/rds/adni_dxs_imputed.rds")
dx.dt |> setkey(PTID, EXAMDATE) |> readr::write_rds(outrds)
rm(outrds)
