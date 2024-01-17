#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
library(Boruta)


## Redo algorithm
reselect_features <- TRUE

## Load data
# baseline data
fpath       <- here("data/rds/triad.rds")
if (file.exists(fpath)) {
  triad.dt  <- read_rds(fpath)
} else {
  here("code/demographics.R") |> source()
}

# Tau and AB accumulation on CEREBRA rois
fpath       <- here("data/pet_biomarkers_cerebra.csv")
if (file.exists(fpath)) {
  cerebra.dt  <- fread(fpath)
} else {
  here("code/parse_csv_data.R") |> source()
}

## Data Cleaning
# Filter
triad.dt    <- triad.dt[!is.na(MOCA_score) &
                        !DX_clean %in% c("Young", "Other", "AD"),
                        #!DX_clean %in% c("Young", "Other"),
                        .(PTID, VISIT, MOCA_score)]

# Cerebra dictionary
# Use the right-side label
dict_roi    <- unique(cerebra.dt[SIDE == "R", .(LABEL_id, LABEL_name)])

# Amyloid
amy_roi.dt  <- cerebra.dt[!is.na(AMYLOID),
                          .(PTID, VISIT, LABEL_name, SIDE,
                            AMYLOID = AMYLOID / VOL)] |>
  dcast(... ~ SIDE, value.var = "AMYLOID")

amy_roi.dt  <- dict_roi[amy_roi.dt, on = "LABEL_name",
                        .(PTID, VISIT,
                          ROI = sprintf("AMY_%03i", LABEL_id),
                          AMYLOID = (L + R) / 2)] |>
  dcast(... ~ ROI, value.var = "AMYLOID")

amy_roi.dt  <- amy_roi.dt[triad.dt, on = .(PTID, VISIT)]
amy_roi.dt  <- amy_roi.dt[complete.cases(amy_roi.dt)]

# Tau
tau_roi.dt  <- cerebra.dt[!is.na(TAU),
                          .(PTID, VISIT, LABEL_name, SIDE,
                            TAU = TAU / VOL)] |>
  dcast(... ~ SIDE, value.var = "TAU")

tau_roi.dt  <- dict_roi[tau_roi.dt, on = "LABEL_name",
                        .(PTID, VISIT,
                          ROI = sprintf("TAU_%03i", LABEL_id),
                          TAU = (L + R) / 2)] |>
  dcast(... ~ ROI, value.var = "TAU")

tau_roi.dt  <- tau_roi.dt[triad.dt, on = .(PTID, VISIT)]
tau_roi.dt  <- tau_roi.dt[complete.cases(tau_roi.dt)]

## Feature selection
# Use the Boruta algorithm to select the most important ROIs related
# Amyloid by Braak Staging
fpath       <- here("data/rds/cerebra_rois_amy_moca.rds")
if (!file.exists(fpath) | reselect_features) {
  set.seed(666)

  cerebra_rois <- stringr::str_subset(names(amy_roi.dt), "AMY")

  # Stage 0 : Lack of tau aggregation
  rois_amy  <- Boruta(amy_roi.dt[, ..cerebra_rois],
                      amy_roi.dt[, MOCA_score]) |>
    #TentativeRoughFix() |>
    #getSelectedAttributes()
    attStats() |>
    as.data.table(keep.rownames = "id")

  write_rds(rois_amy, fpath)
} else {
  rois_amy  <- read_rds(fpath)
}
rm(fpath)

# Tau (all subjects)
fpath       <- here("data/rds/cerebra_rois_tau_moca.rds")
if (!file.exists(fpath) | reselect_features) {
  set.seed(666)

  cerebra_rois <- stringr::str_subset(names(tau_roi.dt), "TAU_\\d{3}")

  rois_tau  <- Boruta(tau_roi.dt[, ..cerebra_rois],
                      tau_roi.dt[, MOCA_score]) |>
    #TentativeRoughFix() |>
    #getSelectedAttributes()
    attStats() |>
    as.data.table(keep.rownames = "id")

  write_rds(rois_tau, fpath)
} else {
  rois_tau  <- read_rds(fpath)
}
rm(fpath)
