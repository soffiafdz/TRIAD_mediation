#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lubridate)
library(ADNIMERGE)


## INPUT
# Raket DTimeline
fpaths      <- here("data", c("adni_raket_cdrsb_mmse.csv",
                              "rds/adni_age_calculated.rds"))

if (!file.exists(fpaths[1])) {
  sprintf("File: %s is required but could not be found.", fpaths[1]) |> stop()
}

raket.dt    <- fread(fpaths[1])

# AGE (for merging)
# Calculated age using the date of MRI sessions
if (file.exists(fpaths[2])) {
  age.dt    <- readRDS(fpaths[2])
} else {
  here("code/calculate_age_adni.R") |> source()
}

rm(fpaths)

# Adnimerge
data(adnimerge)
setDT(adnimerge)

## MERGE
# PTID & EXAMDATE of my clean data
# Include date to avoid future work
id_date.dt  <- adnimerge[, .(RID, Month), .(PTID, VISCODE)
                         ][age.dt, on = .(PTID, VISCODE)
                         ][!is.na(RID),
                         .(PTID, EXAMDATE, AGE, RID, Month = factor(Month))]

# PTID
id.dt       <- id_date.dt[, .(PTID, RID)] |> unique()
raket.dt    <- id.dt[raket.dt[, .(RID, Month = factor(Month), EDT)],
                     on = "RID"][, -"RID"]

# EXAMDATE
raket.dt    <- id_date.dt[raket.dt, on = .(PTID, Month)][, -"RID"]
rm(id.dt, id_date.dt)

# IMPUTE examdate
raket.dt[, EXAMDATE := ymd(EXAMDATE)]

# Filter usable rows (at least one date)
raket.dt    <- raket.dt[raket.dt[!is.na(EXAMDATE), unique(PTID)],
                        on = "PTID"] |> unique()

# Fill out NAs with closest date
setorder(raket.dt, PTID, Month)
raket.dt[, MISS_date := is.na(EXAMDATE)]
raket.dt[, Month := as.integer(as.character(Month))]
raket.dt    <- raket.dt[(!MISS_date), .(PTID, Month, .SD), .SDcols = 2:4
                        ][raket.dt, on = .(PTID, Month), roll = "nearest"] |>
setnames(3:5, c("DATE.src", "AGE.src", "Month.src"))

## Impute data (EXAMDATE & AGE)
raket.dt[(MISS_date), DIFF := Month.src - Month]
raket.dt[(MISS_date),
         EXAMDATE := ymd(round_date(DATE.src - dmonths(DIFF), unit = "day"))]
raket.dt[(MISS_date), AGE := AGE.src - DIFF / 12] # Convert months to years
raket.dt    <- raket.dt[, .(MONTH = Month, AGE, EDT), .(PTID, EXAMDATE)]

### OUTPUT
outrds      <- here("data/rds/adni_raket_edt.rds")
readr::write_rds(raket.dt, outrds)
rm(outrds)
