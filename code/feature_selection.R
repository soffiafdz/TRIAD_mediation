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
                        c("pet_cerebra.rds", "raket_eds.rds", "incl_subs.rds"))

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
if (!file.exists(fpaths[2]))   here("code/parse_csv_data.R") |> source()
if (!file.exists(fpaths[3]))   here("code/demographics.R") |> source()

raket.dt        <- read_rds(fpaths[2]) |> setkey(PTID, VISIT)
all_subs.dt     <- read_rds(fpaths[3]) |> setkey(PTID, VISIT)

amy.dt          <- raket.dt[, .(PTID, VISIT, RAKET_edt, RAKET_group)
                            ][amy.dt][all_subs.dt]

tau.dt          <- raket.dt[, .(PTID, VISIT, RAKET_edt, RAKET_group)
                            ][tau.dt][all_subs.dt]

rm(fpaths, raket.dt, pet.dt, all_subs.dt)

## Feature Selection and Clustering
# Amyloid
fpath           <- here("data/rds/cerebra_rois_raket_amy_all.rds")
if (!file.exists(fpath)) reselect_rois <- TRUE

if (reselect_rois) {
  # Use the Boruta algorithm to select the most important ROIs related
  set.seed(666)

  cerebra_rois  <- str_subset(names(amy.dt), "AMY")

  rois_amy.dt   <-
    Boruta(amy.dt[RAKET_group != "Healthy", ..cerebra_rois],
           amy.dt[RAKET_group != "Healthy", RAKET_edt]) |>
      #TentativeRoughFix() |>
      #getSelectedAttributes()
      attStats() |>
      as.data.table(keep.rownames = "id") |>
      {\(x) x[!decision == "Rejected"]} ()

  rois_amy.dt[, `:=`(id = str_remove(id, "^.{4}"),
                     decision = factor(decision,
                                       levels = c("Confirmed", "Tentative")))]

  write_rds(rois_amy.dt, fpath)
} else {
  rois_amy.dt   <- read_rds(fpath)
}

## Tau
fpath           <- here("data/rds/cerebra_rois_raket_tau_all.rds")
if (!file.exists(fpath)) reselect_rois <- TRUE

if (reselect_rois) {
  # Use the Boruta algorithm to select the most important ROIs related
  set.seed(666)

  cerebra_rois  <- str_subset(names(tau.dt), "TAU")

  rois_tau.dt <-
    Boruta(tau.dt[RAKET_group != "Healthy", ..cerebra_rois],
           tau.dt[RAKET_group != "Healthy", RAKET_edt]) |>
    #TentativeRoughFix() |>
    #getSelectedAttributes()
    attStats() |>
    as.data.table(keep.rownames = "id") |>
    {\(x) x[!decision == "Rejected"]} ()

  rois_tau.dt[, `:=`(id = str_remove(id, "^.{4}"),
                     decision = factor(decision,
                                       levels = c("Confirmed", "Tentative")))]

  write_rds(rois_tau.dt, fpath)
} else {
  rois_tau.dt  <- read_rds(fpath)
}

## Merge ROIS
fname           <- here("data/rds/cerebra_rois_raket_all.rds")
#if (reselect_rois | !file.exists(fname)) {
if (TRUE) {
  rois_amy.dt[, suvr := "amy"]
  rois_tau.dt[, suvr := "tau"]
  rois.dt       <- rbindlist(list(rois_amy.dt, rois_tau.dt))

  #setkey(rois.dt, id)

  #rois_both.dt  <- rois.dt[, .N, .(id, group) ][N == 2, -"N"]

  #rois.dt       <- rbindlist(list(rois_both.dt,
                                  #rois.dt[!rois_both.dt]),
                             #fill = TRUE)

  #rois.dt[is.na(suvr), suvr := "both"]
  #rm(rois_both.dt)
  write_rds(rois.dt, fname)
} else {
  rois.dt       <- read_rds(fname)
}

## Plot
rois.dt |>
  ggplot(aes(x = id, y = meanImp, shape = decision)) +
    theme_classic(base_size = 12) +
    geom_errorbar(aes(ymin = minImp, ymax = maxImp),
                  width = .2, position = position_dodge(.5)) +
    geom_point(fill = "white", size = 3) +
    scale_shape_manual(values = 24:25, guide = "none") +
    labs(x = "Cerebra ROIs", y = "Importance", shape = "Decision") +
    coord_flip() +
    facet_wrap(vars(suvr))

here("plots/boruta_rois_raket_all.png") |>
  ggsave(width = 11, height = 11, units = "in", dpi = 600)
