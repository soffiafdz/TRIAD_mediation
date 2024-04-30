#!/usr/bin/env Rscript

library(here)
library(data.table)
library(stringr)
library(readr)
library(Boruta)
library(ggplot2)
library(ggforce)
library(patchwork)
#library(caret)
#library(cluster)


## Redo algorithm
reselect_rois   <- FALSE

## Read/Parse CSV files
fpaths          <- here("data/rds",
                        c("pet_cerebra.rds", "covars.rds",
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
if (!file.exists(fpaths[4]))   here("code/demographics.R") |> source()
if (any(!file.exists(fpaths[2:3])))  here("code/parse_csv_data.R") |> source()

covars.dt   <- read_rds(fpaths[2]) |> setkey(PTID, VISIT)
raket.dt        <- read_rds(fpaths[3]) |> setkey(PTID, VISIT)
all_subs.dt     <- read_rds(fpaths[4]) |> setkey(PTID, VISIT)

amy.dt          <- covars.dt[, .(PTID, VISIT, MMSE = as.numeric(MMSE))
                             ][raket.dt[, .(PTID, VISIT, RAKET_group)]
                             ][amy.dt][all_subs.dt][!is.na(MMSE)]

tau.dt          <- covars.dt[, .(PTID, VISIT, MMSE = as.numeric(MMSE))
                             ][raket.dt[, .(PTID, VISIT, RAKET_group)]
                             ][tau.dt][all_subs.dt][!is.na(MMSE)]

rm(fpaths, covars.dt, raket.dt, pet.dt, all_subs.dt)

## Feature Selection and Clustering
# Amyloid
fpath           <- here("data/rds/cerebra_rois_mmse_amy.rds")
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
             amy.dt[RAKET_group == raket_grps[i], MMSE]) |>
      #TentativeRoughFix() |>
      #getSelectedAttributes()
      attStats() |>
      as.data.table(keep.rownames = "id") |>
      {\(x) x[!decision == "Rejected"]} ()
  }

  rois_amy.dt   <- rbindlist(rois_amy.lst, idcol = "group")
  rois_amy.dt[, `:=`(id = str_remove(id, "^.{4}"),
                     group = factor(group, levels = c("Healthy",
                                                      "Early stages",
                                                      "Late stages")),
                     decision = factor(decision,
                                       levels = c("Confirmed", "Tentative")))]
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
fpath           <- here("data/rds/cerebra_rois_mmse_tau.rds")
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
             tau.dt[RAKET_group == raket_grps[i], MMSE]) |>
      #TentativeRoughFix() |>
      #getSelectedAttributes()
      attStats() |>
      as.data.table(keep.rownames = "id") |>
      {\(x) x[!decision == "Rejected"]} ()
  }

  rois_tau.dt   <- rbindlist(rois_tau.lst, idcol = "group")
  rois_tau.dt[, `:=`(id = str_remove(id, "^.{4}"),
                     group = factor(group, levels = c("Healthy",
                                                      "Early stages",
                                                      "Late stages")),
                     decision = factor(decision,
                                       levels = c("Confirmed", "Tentative")))]
  rm(rois_tau.lst)

  write_rds(rois_tau.dt, fpath)
} else {
  rois_tau.lst  <- read_rds(fpath)
}

## Plots
p1 <-
  rois_amy.dt |>
  ggplot(aes(x = id, y = meanImp, shape = decision)) +
    theme_classic(base_size = 12) +
    geom_errorbar(aes(ymin = minImp, ymax = maxImp),
                  width = .2, position = position_dodge(.5)) +
    geom_point(fill = "white", size = 3) +
    scale_shape_manual(values = 24:25, guide = "none") +
    labs(x = "Cerebra ROIs", y = "Importance", shape = "Decision",
         title = "Selected ROIs: Amyloid") +
    coord_flip() +
    facet_col(vars(group), scales = "free_y", space = "free")

p2 <-
  rois_tau.dt |>
  ggplot(aes(x = id, y = meanImp, shape = decision)) +
    theme_classic(base_size = 12) +
    geom_errorbar(aes(ymin = minImp, ymax = maxImp),
                  width = .2, position = position_dodge(.5)) +
    geom_point(fill = "white", size = 3) +
    scale_shape_manual(values = 24:25) +
    labs(x = "Cerebra ROIs", y = "Importance", shape = "Decision",
         title = "Selected ROIs: Tau",
         caption = "Features selected to predict MMSE using Boruta") +
    coord_flip() +
    facet_col(vars(group), scales = "free_y", space = "free")

pp <- p1 + p2
here("plots/boruta-rois_mmse.png") |>
  ggsave(pp, width = 11, height = 11, units = "in", dpi = 600)

## Find selected ROIs that are present on both Amy and Tau lists
if (reselect_rois) {
  rois.dt       <- rbindlist(list(rois_amy.dt[, .(id, group, suvr = "amy")],
                                  rois_tau.dt[, .(id, group, suvr = "tau")]))

  setkey(rois.dt, id, group)

  rois_both.dt  <- rois.dt[, .N, .(id, group) ][N == 2, -"N"]

  rois.dt       <- rbindlist(list(rois_both.dt,
                                  rois.dt[!rois_both.dt]),
                             fill = TRUE)

  rois.dt[is.na(suvr), suvr := "both"]
  write_rds(rois.dt, here("data/rds/cerebra_rois_mmse.rds"))
  rm(rois_both.dt)
}
