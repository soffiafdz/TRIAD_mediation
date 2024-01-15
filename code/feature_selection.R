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
# Convert Braak Staging to float
triad.dt[, TAU_braak_stage := as.numeric(TAU_braak_stage)]

# Filter
triad.dt    <- triad.dt[!is.na(MOCA_score) &
                        !is.na(TAU_braak_group) &
                        !DX_clean %in% c("Young", "Other"),
                        .(PTID, VISIT, TAU_braak_group, MOCA_score)]

# Cerebra dictionary
dict_roi    <- unique(cerebra.dt[, .(LABEL_id, LABEL_name, SIDE)])

# Amyloid
amy_roi.dt  <- cerebra.dt[!is.na(AMYLOID),
                          .(PTID, VISIT, AMYLOID,
                            ROI = sprintf("AMY_%03i", LABEL_id))] |>
  dcast(... ~ ROI, value.var = "AMYLOID")

amy_roi.dt  <- amy_roi.dt[triad.dt, on = .(PTID, VISIT)]
amy_roi.dt  <- amy_roi.dt[complete.cases(amy_roi.dt)]

# Tau
tau_roi.dt  <- cerebra.dt[!is.na(TAU),
                          .(PTID, VISIT, TAU,
                            ROI = sprintf("TAU_%03i", LABEL_id))] |>
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
  rois_amy_0  <- Boruta(amy_roi.dt[TAU_braak_group == "0", ..cerebra_rois],
                        amy_roi.dt[TAU_braak_group == "0", MOCA_score]) |>
    TentativeRoughFix() |>
    #getSelectedAttributes()
    attStats() |>
    as.data.table(keep.rownames = "id")

  # Stage 1 & 2 : Initial stages
  rois_amy_12 <- Boruta(amy_roi.dt[TAU_braak_group == "1 & 2", ..cerebra_rois],
                        amy_roi.dt[TAU_braak_group == "1 & 2", MOCA_score]) |>
    TentativeRoughFix() |>
    #getSelectedAttributes()
    attStats() |>
    as.data.table(keep.rownames = "id")

  # Stage 3 & 4 : Moderate accumulation
  rois_amy_34 <- Boruta(amy_roi.dt[TAU_braak_group == "3 & 4", ..cerebra_rois],
                        amy_roi.dt[TAU_braak_group == "3 & 4", MOCA_score]) |>
    TentativeRoughFix() |>
    #getSelectedAttributes()
    attStats() |>
    as.data.table(keep.rownames = "id")

  # Stage 5 & 6 : Severe accumulation
  rois_amy_56 <- Boruta(amy_roi.dt[TAU_braak_group == "5 & 6", ..cerebra_rois],
                        amy_roi.dt[TAU_braak_group == "5 & 6", MOCA_score]) |>
    TentativeRoughFix() |>
    #getSelectedAttributes()
    attStats() |>
    as.data.table(keep.rownames = "id")

  rois_amy  <- list(rois_amy_0, rois_amy_12, rois_amy_34, rois_amy_56)
  names(rois_amy) <- paste0("Braak_stage_", c("0", "1&2", "3&4", "5&6"))
  write_rds(rois_amy, fpath)
  rm(rois_amy_0, rois_amy_12, rois_amy_34, rois_amy_56)
} else {
  rois_amy  <- read_rds(fpath)
}
rm(fpath)

# Tau (all subjects)
fpath       <- here("data/rds/cerebra_rois_tau_moca.rds")
if (!file.exists(fpath) | reselect_features) {
  set.seed(666)

  cerebra_rois <- stringr::str_subset(names(tau_roi.dt), "TAU_\\d{3}")

  #rois_tau  <- Boruta(tau_roi.dt[, ..cerebra_rois],
                      #tau_roi.dt[, MOCA_score]) |>
    #TentativeRoughFix() |>
    ##getSelectedAttributes()
    #attStats() |>
    #as.data.table(keep.rownames = "id")

  #write_rds(rois_tau, fpath)

  # Stage 0 : Lack of tau aggregation
  rois_tau_0  <- Boruta(tau_roi.dt[TAU_braak_group == "0", ..cerebra_rois],
                        tau_roi.dt[TAU_braak_group == "0", MOCA_score]) |>
    TentativeRoughFix() |>
    #getSelectedAttributes()
    attStats() |>
    as.data.table(keep.rownames = "id")

  # Stage 1 & 2 : Initial stages
  rois_tau_12 <- Boruta(tau_roi.dt[TAU_braak_group == "1 & 2", ..cerebra_rois],
                        tau_roi.dt[TAU_braak_group == "1 & 2", MOCA_score]) |>
    TentativeRoughFix() |>
    #getSelectedAttributes()
    attStats() |>
    as.data.table(keep.rownames = "id")

  # Stage 3 & 4 : Moderate accumulation
  rois_tau_34 <- Boruta(tau_roi.dt[TAU_braak_group == "3 & 4", ..cerebra_rois],
                        tau_roi.dt[TAU_braak_group == "3 & 4", MOCA_score]) |>
    TentativeRoughFix() |>
    #getSelectedAttributes()
    attStats() |>
    as.data.table(keep.rownames = "id")

  # Stage 5 & 6 : Severe accumulation
  rois_tau_56 <- Boruta(tau_roi.dt[TAU_braak_group == "5 & 6", ..cerebra_rois],
                        tau_roi.dt[TAU_braak_group == "5 & 6", MOCA_score]) |>
    TentativeRoughFix() |>
    #getSelectedAttributes()
    attStats() |>
    as.data.table(keep.rownames = "id")

  rois_tau  <- list(rois_tau_0, rois_tau_12, rois_tau_34, rois_tau_56)
  names(rois_tau) <- paste0("Braak_stage_", c("0", "1&2", "3&4", "5&6"))
  write_rds(rois_tau, fpath)
  rm(rois_tau_0, rois_tau_12, rois_tau_34, rois_tau_56)
} else {
  rois_tau  <- read_rds(fpath)
}
rm(fpath)
