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
recluster       <- FALSE

## Read/Parse CSV files
fpaths          <- here("data/rds",
                        c("pet_cerebra.rds", #"covars.rds",
                          "raket_eds.rds", "incl_subs.dt"))

## Preprocess PET data
if (!file.exists(fpaths[1]))   here("code/parse_pet.R") |> source()
pet.dt          <- read_rds(fpaths[1]) |> setkey(PTID, VISIT)

# Amyloid (NAV)
amy.dt          <- pet.dt[, .(PTID, VISIT,
                              ROI = paste("AMY", LABEL_id, sep = "_"),
                              SUVR_norm_log = log(SUVR_nav_norm))] |>
  dcast(... ~ ROI, value.var = "SUVR_norm_log")


## Tau (MK)
tau.dt          <- pet.dt[, .(PTID, VISIT,
                              ROI = paste("TAU", LABEL_id, sep = "_"),
                              SUVR_norm_log = log(SUVR_mk_norm))] |>
  dcast(... ~ ROI, value.var = "SUVR_norm_log")

## Covariates
if (!file.exists(fpaths[3]))   here("code/demographics.R") |> source()
#if (any(!file.exists(fpaths[2:3])))  here("code/parse_csv_data.R") |> source()

#covars.dt   <- read_rds(fpaths[2]) |> setkey(PTID, VISIT)
raket.dt        <- read_rds(fpaths[2]) |> setkey(PTID, VISIT)
all_subs.dt     <- read_rds(fpaths[3]) |> setkey(PTID, VISIT)

amy.dt          <- raket.dt[, .(PTID, VISIT, RAKET_edt, RAKET_group)
                            ][amy.dt][all_subs.dt]

tau.dt          <- raket.dt[, .(PTID, VISIT, RAKET_edt, RAKET_group)
                            ][tau.dt][all_subs.dt]

#rm(fpaths, covars.dt, raket.dt, pet.dt, all_subs.dt)

## Feature Selection and Clustering
# Amyloid
fpath           <- here("data/rds/cerebra_rois_amy_raket.rds")
if (!file.exists(fpath)) reselect_rois <- TRUE

if (reselect_rois) {
  # Use the Boruta algorithm to select the most important ROIs related
  set.seed(666)

  cerebra_rois  <- str_subset(names(amy.dt), "AMY")

  raket_grps    <- amy.dt[, levels(RAKET_group)]
  rois_amy.lst  <- vector("list", 3)
  names(rois_amy.lst)   <- raket_grps

  for (i in seq_along(raket_grps)) {
    rois_amy.lst[[i]]   <-
      Boruta(amy.dt[RAKET_group == raket_grps[i], ..cerebra_rois],
             amy.dt[RAKET_group == raket_grps[i], RAKET_edt]) |>
      #TentativeRoughFix() |>
      #getSelectedAttributes()
      attStats() |>
      as.data.table(keep.rownames = "id") |>
      {\(x) x[!decision == "Rejected"]} ()
  }

  rois_amy.dt   <- rbindlist(rois_amy.lst, idcol = "group")
  rm(rois_amy.lst)

  write_rds(rois_amy.dt, fpath)
} else {
  rois_amy.dt   <- read_rds(fpath)
}

##if (FALSE) {
##if (reselect_rois | recluster) {
  ## Clustering K-means
  #K             <- 5

  #cerebra_rois  <- str_subset(names(amy.dt), "AMY")
  #amy_scaled.dt <- amy.dt[, lapply(.SD, scale), .SDcols = cerebra_rois]
  #names(amy_scaled.dt)  <- str_remove(names(amy_scaled.dt), ".{3}$")

  ## Only selected ROIs from Boruta
  #rois_amy      <- rois_amy.dt[decision == "Confirmed", id]
  #amy_rois_scl.dt       <- amy_scaled.dt[, ..rois_amy]

  #set.seed(42)
  #clusters_amy  <- kmeans(transpose(amy_rois_scl.dt), centers = K)
  #rois_amy.dt[decision == "Confirmed", cluster := clusters_amy$cluster]
#}

#if (reselect_rois | recluster) write_rds(rois_amy.dt, fpath)
#rm(fpath)

## Tau
fpath           <- here("data/rds/cerebra_rois_tau_raket.rds")
if (!file.exists(fpath)) reselect_rois <- TRUE

if (reselect_rois) {
  # Use the Boruta algorithm to select the most important ROIs related
  set.seed(666)

  cerebra_rois  <- str_subset(names(tau.dt), "TAU")

  raket_grps    <- tau.dt[, levels(RAKET_group)]
  rois_tau.lst  <- vector("list", 3)
  names(rois_tau.lst)   <- raket_grps

  for (i in seq_along(raket_grps)) {
    rois_tau.lst[[i]]   <-
      Boruta(tau.dt[RAKET_group == raket_grps[i], ..cerebra_rois],
             tau.dt[RAKET_group == raket_grps[i], RAKET_edt]) |>
      #TentativeRoughFix() |>
      #getSelectedAttributes()
      attStats() |>
      as.data.table(keep.rownames = "id") |>
      {\(x) x[!decision == "Rejected"]} ()
  }

  rois_tau.dt   <- rbindlist(rois_tau.lst, idcol = "group")
  rm(rois_tau.lst)

  write_rds(rois_tau.dt, fpath)
} else {
  rois_tau.lst  <- read_rds(fpath)
}

##if (FALSE) {
##if (reselect_rois | recluster) {
  ## Clustering K-means
  #K             <- 5

  #cerebra_rois  <- str_subset(names(tau.dt), "TAU")
  #tau_scaled.dt <- tau.dt[, lapply(.SD, scale), .SDcols = cerebra_rois]
  #names(tau_scaled.dt)  <- str_remove(names(tau_scaled.dt), ".{3}$")

  ## Only selected ROIs from Boruta
  #rois_tau      <- rois_tau.dt[decision == "Confirmed", id]
  #tau_rois_scl.dt       <- tau_scaled.dt[, ..rois_tau]

  #set.seed(42)
  #clusters_tau  <- kmeans(transpose(tau_rois_scl.dt), centers = K)
  #rois_tau.dt[decision == "Confirmed", cluster := clusters_tau$cluster]
#}

#if (reselect_rois | recluster) write_rds(rois_tau.dt, fpath)
#rm(fpath)

## Find selected ROIs that are present on both Amy and Tau lists
#if (reselect_rois) {
  ## In ANY list
  #rois_amy_tau_a  <- c(rois_amy.dt[decision == "Confirmed",
                                   #str_extract(id, "\\d{3}")],
                       #rois_tau.dt[decision == "Confirmed",
                                   #str_extract(id, "\\d{3}")]) |>
  #unique() |> as.numeric()

  ## In BOTH lists
  #rois_amy_tau_b  <- rois_amy.dt[decision == "Confirmed",
                                 #.(str_extract(id, "\\d{3}"))
                                 #][rois_tau.dt[decision == "Confirmed",
                                               #.(str_extract(id, "\\d{3}"))],
                                 #on = "V1", nomatch = 0, as.numeric(V1)]

  #rois_amy_tau.dt <-
    #dict_roi[LABEL_id %in% rois_amy_tau_a][order(LABEL_name)]

  #rois_amy_tau.dt[, LIST := "ANY"]
  #rois_amy_tau.dt[LABEL_id %in% rois_amy_tau_b, LIST := "BOTH"]

  #write_rds(rois_amy_tau.dt, here("data/rds/cerebra_rois_amy_tau_moca.rds"))
  #rm(rois_amy_tau_a, rois_amy_tau_b)
#}
