#!/usr/bin/env Rscript

library(here)
library(data.table)
library(stringr)
library(readr)
library(Boruta)
library(caret)
library(cluster)


## Redo algorithm
reselect_rois   <- TRUE
recluster       <- TRUE

## Load data
# baseline data
fpath           <- here("data/rds/triad.rds")
if (file.exists(fpath)) {
  triad.dt      <- read_rds(fpath)
} else {
  here("code/demographics.R") |> source()
}

# Tau and AB accumulation on CEREBRA rois
fpath           <- here("data/pet_biomarkers_cerebra.csv")
if (file.exists(fpath)) {
  cerebra.dt    <- fread(fpath)
} else {
  here("code/parse_csv_data.R") |> source()
}

## Data Cleaning
# Filter
triad.dt        <- triad.dt[!is.na(MOCA_score) &
                            !DX_clean %in% c("Young", "Other", "AD"),
                            #!DX_clean %in% c("Young", "Other"),
                            .(PTID, VISIT, MOCA_score)]

# Cerebra dictionary
# Use the right-side label
dict_roi        <- unique(cerebra.dt[, .(LABEL_id, LABEL_name, SIDE)])

# Amyloid
amy.dt          <- cerebra.dt[!is.na(AMYLOID_norm),
                                  .(PTID, VISIT,
                                    ROI = sprintf("AMY_%03i", LABEL_id),
                                    #AMYLOID = log(AMYLOID_norm / VOL))] |>
                                    AMYLOID = AMYLOID_norm / VOL)] |>
  dcast(... ~ ROI, value.var = "AMYLOID")

amy.dt          <- amy.dt[triad.dt, on = .(PTID, VISIT)]
amy.dt          <- amy.dt[complete.cases(amy.dt)]

# Tau
tau.dt          <- cerebra.dt[!is.na(TAU_norm),
                                  .(PTID, VISIT,
                                    ROI = sprintf("TAU_%03i", LABEL_id),
                                    #TAU = log(TAU_norm / VOL))] |>
                                    TAU = TAU_norm / VOL)] |>
  dcast(... ~ ROI, value.var = "TAU")

tau.dt          <- tau.dt[triad.dt, on = .(PTID, VISIT)]
tau.dt          <- tau.dt[complete.cases(tau.dt)]


## Feature Selection and Clustering
fpath           <- here("data/rds/cerebra_rois_amy_moca.rds")
if (!file.exists(fpath)) reselect_rois <- TRUE

if (reselect_rois) {
  # Use the Boruta algorithm to select the most important ROIs related
  set.seed(666)

  cerebra_rois  <- str_subset(names(amy.dt), "AMY")

  rois_amy.dt   <- Boruta(amy.dt[, ..cerebra_rois],
                      amy.dt[, MOCA_score]) |>
    #TentativeRoughFix() |>
    #getSelectedAttributes()
    attStats() |>
    as.data.table(keep.rownames = "id")
} else {
  rois_amy.dt   <- read_rds(fpath)
}

if (reselect_rois | recluster) {
  # Clustering K-means
  K               <- 5

  amy_cols        <- names(amy.dt)[startsWith(names(amy.dt), "AMY")]
  amy_scaled.dt   <- amy.dt[, lapply(.SD, scale), .SDcols = amy_cols]
  names(amy_scaled.dt) <- str_remove(names(amy_scaled.dt), ".{3}$")

  # Only selected ROIs from Boruta
  rois_amy        <- rois_amy.dt[decision == "Confirmed", id]
  amy_rois_scl.dt <- amy_scaled.dt[, ..rois_amy]

  set.seed(42)
  clusters_amy    <- kmeans(transpose(amy_rois_scl.dt), centers = K)
  rois_amy.dt[decision == "Confirmed", cluster := clusters_amy$cluster]
}

if (reselect_rois | recluster) write_rds(rois_amy.dt, fpath)
rm(fpath)

# Tau (all subjects)
fpath           <- here("data/rds/cerebra_rois_tau_moca.rds")
if (!file.exists(fpath)) reselect_rois <- TRUE

if (reselect_rois) {
  # Boruta algorithm
  set.seed(666)

  cerebra_rois  <- str_subset(names(tau.dt), "TAU_\\d{3}")

  rois_tau.dt   <- Boruta(tau.dt[, ..cerebra_rois], tau.dt[, MOCA_score]) |>
    #TentativeRoughFix() |>
    #getSelectedAttributes()
    attStats() |>
    as.data.table(keep.rownames = "id")

  write_rds(rois_tau.dt, fpath)
} else {
  rois_tau.dt   <- read_rds(fpath)
}

if (reselect_rois | recluster) {
  # Clustering K-means
  K               <- 5

  tau_cols        <- names(tau.dt)[startsWith(names(tau.dt), "TAU")]
  tau_scaled.dt   <- tau.dt[, lapply(.SD, scale), .SDcols = tau_cols]
  names(tau_scaled.dt) <- str_remove(names(tau_scaled.dt), ".{3}$")

  # Only selected ROIs from Boruta
  rois_tau        <- rois_tau.dt[decision == "Confirmed", id]
  tau_rois_scl.dt <- tau_scaled.dt[, ..rois_tau]

  set.seed(42)
  clusters_tau    <- kmeans(transpose(tau_rois_scl.dt), centers = K)
  rois_tau.dt[decision == "Confirmed", cluster := clusters_tau$cluster]
}

if (reselect_rois | recluster) write_rds(rois_tau.dt, fpath)
rm(fpath)

# Find selected ROIs that are present on both Amy and Tau lists
if (reselect_rois) {
  rois_amy_tau    <-
    dict_roi[rois_amy.dt[decision == "Confirmed",
                         .(id = as.numeric(str_extract(id, "\\d{3}")))
             ][rois_tau.dt[decision == "Confirmed",
                           .(id = as.numeric(str_extract(id, "\\d{3}")))
             ], on = "id", nomatch = 0] , on = .(LABEL_id = id)
             ][order(LABEL_name)]

  write_rds(rois_amy_tau, here("data/rds/cerebra_rois_amy_tau_moca.rds"))
}
