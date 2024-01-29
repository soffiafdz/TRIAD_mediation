#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
library(stringr)
library(progress)
library(lavaan)
library(lavaanPlot)
library(ggplot2)
library(patchwork)
#library(semTable)

## Refit models
refit_mods  <- FALSE

## Print plots ### Needs to be done outside renv
print_plots <- TRUE

## INPUT
# Load baseline data
fpath       <- here("data/rds/triad.rds")
if (file.exists(fpath)) {
  triad.dt  <- read_rds(fpath)
} else {
  here("code/demographics.R") |> source()
}
rm(fpath)

# Load Cerebra data
# Biomarkers
fpath       <- here("data/pet_biomarkers_cerebra.csv")
if (file.exists(fpath)) {
  cerebra.dt      <- fread(fpath)
} else {
  here("code/parse_csv_data.R") |> source()
}

# ROIs
fpaths      <- here("data/rds",
                    sprintf("cerebra_rois_%s_moca.rds",
                            c("amy", "tau", "amy_tau")))
if (any(!file.exists(fpaths))) {
  here("code/feature_selection.R") |> source()
} else {
  rois_amy.dt     <- read_rds(fpaths[1])
  rois_tau.dt     <- read_rds(fpaths[2])
  rois_amy_tau.dt <- read_rds(fpaths[3])
}
rm(fpaths)

## Data cleaning
# Convert Sex to dummy variable
triad.dt[, SEX_n := as.numeric(SEX) - 1]

# Remove youth and dementias
triad.dt    <- triad.dt[!DX_clean %in% c("Young", "Other", "AD")]

# 1 - HVR (average for both sides)
#triad.dt[, `:=`(HVR_lr = 1 - HVR_l, HVR_rr = 1 - HVR_r)]
triad.dt[, `:=`(HVR_mean_inv = 1 - (HVR_l + HVR_r) / 2)]

triad.dt    <- triad.dt[, .(PTID, VISIT, DX, SEX_n, AGE_scan, EDUC,#APOE_n,
                            HVR_mean_inv, MOCA_score)]

# Merge triad and cerebra data
# Cerebra dictionary
dict_roi    <- unique(cerebra.dt[, .(LABEL_id, LABEL_name, SIDE)])

# Selected features for both Amyloid and Tau
# Weighted average of all ROIs
#amy.dt      <- cerebra.dt[LABEL_id %in% rois_amy_tau.dt$LABEL_id
amy.dt      <- cerebra.dt[LABEL_id %in% rois_amy_tau.dt[LIST == "BOTH",
                                                        LABEL_id]
                          ][!is.na(AMYLOID_norm),
                          .(AMYLOID = weighted.mean(AMYLOID_norm, VOL)),
                          .(PTID, VISIT)]

#tau.dt      <- cerebra.dt[LABEL_id %in% rois_amy_tau.dt$LABEL_id
tau.dt      <- cerebra.dt[LABEL_id %in% rois_amy_tau.dt[LIST == "BOTH",
                                                        LABEL_id]
                          ][!is.na(AMYLOID_norm),
                          .(TAU = weighted.mean(TAU_norm, VOL)),
                          .(PTID, VISIT)]

triad_w.dt  <- tau.dt[amy.dt, on = .(PTID, VISIT)
                      ][triad.dt, on = .(PTID, VISIT)
                      ][!is.na(MOCA_score) & !is.na(AMYLOID) & !is.na(TAU)]

# ROIs normalized by Volume
amy.dt      <- cerebra.dt[LABEL_id %in% rois_amy_tau.dt$LABEL_id
                          ][!is.na(AMYLOID_norm),
                          .(PTID, VISIT, LABEL_id,
                            AMYLOID = AMYLOID_norm / VOL * 1000)]

tau.dt      <- cerebra.dt[LABEL_id %in% rois_amy_tau.dt$LABEL_id
                          ][!is.na(AMYLOID_norm),
                          .(PTID, VISIT, LABEL_id,
                            TAU = TAU_norm / VOL * 1000)]

triad.dt    <- tau.dt[amy.dt, on = .(PTID, VISIT, LABEL_id)
                      ][triad.dt, on = .(PTID, VISIT)
                      ][!is.na(MOCA_score) & !is.na(AMYLOID) & !is.na(TAU)]
rm(amy.dt, tau.dt)

## Model definitions
hvr.mod     <-
  str_glue("
# Regressions
AMYLOID ~ SEX_n + AGE_scan
TAU ~ a * AMYLOID + SEX_n + AGE_scan
HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan
# HVR mediation
# Direct effect
dAMY  := b
# Indirect effect
iTAU  := a * c
# Total effect
Total := dAMY + iTAU
# Proportion analysis
propAMY := dAMY / Total
propTAU := iTAU / Total
")

moca.mod    <-
  str_glue("
# Regressions
AMYLOID ~ SEX_n + AGE_scan
TAU ~ a * AMYLOID + SEX_n + AGE_scan
HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan
MOCA_score ~ d * AMYLOID + e * TAU + f * HVR_mean_inv + SEX_n + AGE_scan + EDUC
# Direct effect
dAMY := d
# Indirect effects
iTAU := a * e
iHVR := (b + (a * c)) * f
# Total effects
iTotal    := iTAU + iHVR
Total := iTotal + dAMY
# Proportion analysis
propAMY := dAMY / Total
propTAU := iTAU / Total
propHVR := iHVR / Total
propInd := iTotal / Total
")


## Fit models
# General models with weighted average of all selected ROIs
fnames          <- here("data/rds", sprintf("mediation_%s_w_bs.rds",
                                            c("hvr", "moca")))

if (any(!file.exists(fnames)) | refit_mods) {
  mods.lst      <- vector("list", 2)
  for (i in seq_along(mods.lst)) {
    mods.lst[[i]]  <-
      sem(c(hvr.mod, moca.mod)[i],
          data = triad_w.dt,
          estimator = "ML",
          se = "bootstrap",
          bootstrap = 1000)

    write_rds(mods.lst[[i]], fnames[i])
  }

  mod_hvr_w.fit         <- mods.lst[[1]]
  mod_moca_w.fit        <- mods.lst[[2]]
  rm(mods.lst)
} else {
  mod_hvr_w.fit         <- read_rds(fnames[1])
  mod_moca_w.fit        <- read_rds(fnames[2])
}
rm(fnames)

# Individual models by ROI
fnames          <- here("data/rds",
                        sprintf("mediation_%s_bs.rds", c("hvr", "moca")))

if (any(!file.exists(fnames)) | refit_mods) {
  mod_hvr.fits  <- mod_moca.fits <- vector("list", rois_amy_tau.dt[, .N])
  mods.lst      <- list(mod_hvr.fits, mod_moca.fits)

  pb <- progress_bar$new(format = "Models | :what [:bar] :current/:total",
                         total = length(mods.lst) * rois_amy_tau.dt[, .N],
                         clear = FALSE, width = 75)

  for (i in seq_along(mods.lst)) {
    names(mods.lst[[i]]) <- rois_amy_tau.dt[, paste(LABEL_name, SIDE,
                                                    sep = "_")]

    for (j in seq_along(rois_amy_tau.dt$LABEL_id)) {
      pb$tick(tokens = list(what = sprintf("%s : %s",
                                           c("HVR", "MoCA")[i],
                                           rois_amy_tau.dt[j, LABEL_name])))

      mods.lst[[i]][[j]]  <-
        sem(c(hvr.mod, moca.mod)[i],
            data = triad.dt[LABEL_id == rois_amy_tau.dt[j, LABEL_id]],
            estimator = "ML",
            se = "bootstrap",
            bootstrap = 1000)
    }

    write_rds(mods.lst[[i]], fnames[i])
  }

  mod_hvr.fits          <- mods.lst[[1]]
  mod_moca.fits         <- mods.lst[[2]]
  rm(mods.lst)
} else {
  mod_hvr.fits          <- read_rds(fnames[1])
  mod_moca.fits         <- read_rds(fnames[2])
}
rm(fnames)

## Extract Fit measures and paramater estimates
# Names for data cleaning
msrs_names      <- mod_hvr.fits[[1]] |> fitMeasures() |> names()
rois_names      <- rois_amy_tau.dt[, paste(LABEL_name, SIDE, sep = "_")]

# HVR model
mod_hvr.msrs    <- mod_hvr.fits |>
  lapply(fitMeasures) |>
  lapply(as.data.table) |>
  lapply(transpose) |>
  rbindlist()

setnames(mod_hvr.msrs, msrs_names)
mod_hvr.msrs[, ROI := rois_names]

mod_hvr.est     <- mod_hvr.fits |>
  lapply(standardizedSolution) |>
  lapply(as.data.table)

mod_hvr.est     <- rois_names |>
  lapply(function (name) mod_hvr.est[[name]][, ROI := name]) |>
  rbindlist()

mod_hvr.est     <- mod_hvr.est[op != "~~"]

mod_moca.msrs   <- mod_moca.fits |>
  lapply(fitMeasures) |>
  lapply(as.data.table) |>
  lapply(transpose) |>
  rbindlist()

setnames(mod_moca.msrs, msrs_names)
mod_hvr.msrs[, ROI := rois_names]

mod_moca.est    <- mod_moca.fits |>
  lapply(standardizedSolution) |>
  lapply(as.data.table)

mod_moca.est    <- rois_names |>
  lapply(function (name) mod_moca.est[[name]][, ROI := name]) |>
  rbindlist()

mod_moca.est    <- mod_moca.est[op != "~~"]

## Plot standardized estimates
# Use only ROIs in BOTH AMY and TAU important lists


# HVR model
DT <- mod_hvr.est[lhs %in% c("dAMY", "iTAU") &
                  ROI %in% rois_amy_tau.dt[LIST == "BOTH",
                                           paste(LABEL_name, SIDE,
                                                 sep = "_")]]
DT[, label := factor(label, levels = c("dAMY", "iTAU"),
                     labels = c("Direct: Amyloid", "Indirect: Tau"))]

ordered_rois  <- DT[lhs == "iTAU"][order(est.std), ROI]

p1 <- DT |>
  ggplot(aes(x = ROI, y = est.std)) +
  theme_classic(base_size = 12) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper),
                width = 0.2, position = position_dodge(0.5)) +
  geom_point(shape = 21, fill = "white", size = 1.5) +
  scale_x_discrete(limits = ordered_rois) +
  labs(x = "CerebrA ROIs", y = "Standardized estimates",
       title = "Direct & indirect effects on HC atrophy",
       caption = "Bootstrap CIs: 1000 resamples") +
  coord_flip() +
  facet_wrap(vars(label))

#here("plots/mediation_hvr_std-estimates.png") |>
  #ggsave(p, width = 6, height = 5, units = "in", dpi = 600)

# MoCA model
#mod_moca.est[lhs == "dAMY", label := "Direct: Amyloid"]
#mod_moca.est[lhs == "iTAU", label := "Indirect: Tau"]
#mod_moca.est[lhs == "iHVR", label := "Indirect: HVR"]

#DT <- mod_moca.est[lhs %in% c("dAMY", "iTAU", "iHVR", "Total")]
#DT[, label := factor(label,
                     #levels = c("Total", "dAMY", "iTAU", "iHVR"),
                     #labels = c("Total", "Direct: Amyloid",
                                #"Indirect: Tau", "Indirect: HVR"))]

DT <- mod_moca.est[lhs %in% c("dAMY", "iTAU") &
                  ROI %in% rois_amy_tau.dt[LIST == "BOTH",
                                           paste(LABEL_name, SIDE,
                                                 sep = "_")]]
DT[, label := factor(label, levels = c("dAMY", "iTAU"),
                     labels = c("Direct: Amyloid", "Indirect: Tau"))]

ordered_rois  <- DT[lhs == "iTAU"][order(-est.std), ROI]

p2 <- DT |>
  ggplot(aes(x = ROI, y = est.std)) +
  theme_classic(base_size = 12) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper),
                width = 0.2, position = position_dodge(0.5)) +
  geom_point(shape = 21, fill = "white", size = 1.5) +
  scale_x_discrete(limits = ordered_rois) +
  labs(x = "CerebrA ROIs", y = "Standardized estimates",
       title = "Direct & indirect effects on MoCA scores",
       caption = "Bootstrap CIs: 1000 resamples") +
  coord_flip() +
  facet_wrap(vars(label))

pp <- p1 + p2
here("plots/mediation_std-estimates.png") |>
  ggsave(pp, width = 12, height = 5, units = "in", dpi = 600)


# Path plot of selected ROIs model
fname       <- here("plots/mediation_moca_path.pdf")
if (!file.exists(fname) | print_plots) {
labels      <- c(AGE_scan = "Age", SEX_n = "Sex", EDUC = "Education",
                 HVR_mean_inv = "HC-atrophy", MOCA_score = "MoCA")

  p_plot    <- lavaanPlot2(model = mod_moca_w.fit, labels = labels,
                            graph_options = list(rankdir = "LR"),
                            node_options = list(shape = "box"),
                            edge_options = list(color = "grey"),
                            coef_labels = T, stand = T,
                            stars = "regress")
  embed_plot_pdf(p_plot, fname)
}
rm(fname)
