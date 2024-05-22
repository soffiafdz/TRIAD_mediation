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
fpaths      <- here("data/rds", c("covars.rds",
                                  "raket_eds.rds",
                                  "wmh_vols.rds",
                                  "hcv_hvr_adj-old.rds",
                                  "pet_cerebra.rds",
                                  "incl_subs.rds",
                                  "cerebra_rois_raket.rds"))

if (!file.exists(fpaths[7]))    here("code/feature_selection.R")  |> source()
if (!file.exists(fpaths[6]))    here("code/demographics.R")       |> source()
if (!file.exists(fpaths[5]))    here("code/parse_pet.R")          |> source()
if (!file.exists(fpaths[4]))    here("code/calc_hvr.R")           |> source()
if (any(!file.exists(fpaths)))  here("code/parse_csv_data.R")     |> source()

covars.dt   <- read_rds(fpaths[1]) |> setkey(PTID, VISIT)
raket.dt    <- read_rds(fpaths[2]) |> setkey(PTID, VISIT)
wmh.dt      <- read_rds(fpaths[3]) |> setkey(PTID, VISIT)
vols.dt     <- read_rds(fpaths[4]) |> setkey(PTID, VISIT)
pet.dt      <- read_rds(fpaths[5]) |> setkey(PTID, VISIT)
all_subs.dt <- read_rds(fpaths[6]) |> setkey(PTID, VISIT)
rois.dt     <- read_rds(fpaths[7])

rm(fpaths)

## Data cleaning
DT <-
  covars.dt[, .(PTID, VISIT, AGE, SEX, EDUC, MMSE)
            ][raket.dt[AB_bool == TRUE, .(PTID, VISIT, RAKET_group)]
            ][wmh.dt
            ][vols.dt[, .(PTID, VISIT,
                          HCv_mean      = (HCv_l + HCv_r) / 2,
                          HVR_mean_inv  = 1 - (HVR_l + HVR_r) / 2)]
            ][all_subs.dt]

# Convert Sex to dummy variable
DT[, SEX_n := as.numeric(SEX) - 1]

# rename
setnames(DT, c("HCv_mean", "HVR_mean_inv"), c("HCv", "HVR"))

# Calculate weighted mean of different networks
rois.dt     <- rois.dt[, .(ids = list(id)),
                       .(suvr, group = factor(group, labels = 0:2))]

networks    <- split(rois.dt$ids, list(rois.dt$suvr, rois.dt$group)) |>
                lapply(unlist)

pet_nets.dt <- unique(pet.dt[, .(PTID, VISIT)])

for (i in rev(seq_along(networks))) {
  net.dt    <- pet.dt[LABEL_id %in% networks[[i]],
                      fifelse(startsWith(names(networks[i]), "amy"),
                              weighted.mean(SUVR_nav, VOLUME_nav),
                              weighted.mean(SUVR_mk, VOLUME_mk)),
                      .(PTID, VISIT)]
  setnames(net.dt, "V1", names(networks[i]))
  DT        <- net.dt[DT]
  rm(net.dt)
}

## Model definitions
nets.mods   <-
  list(str_glue("
# Regressions
amy.0 ~ SEX_n + AGE
tau.0 ~ a * amy.0 + SEX_n + AGE
"),
str_glue("
# Regressions
amy.0 ~ SEX_n + AGE
amy.1 ~ SEX_n + AGE
tau.0 ~ a0 * amy.0 + a1 * amy.1 + SEX_n + AGE
tau.1 ~ b0 * amy.0 + b1 * amy.1 + SEX_n + AGE
"),
str_glue("
# Regressions
amy.0 ~ SEX_n + AGE
amy.1 ~ SEX_n + AGE
amy.2 ~ SEX_n + AGE
tau.0 ~ a0 * amy.0 + a1 * amy.1 + a2 * amy.2 + SEX_n + AGE
tau.1 ~ b0 * amy.0 + b1 * amy.1 + b2 * amy.2 + SEX_n + AGE
tau.2 ~ c0 * amy.0 + c1 * amy.1 + c2 * amy.2 + SEX_n + AGE
"))

hvr.mods    <-
  list(str_glue("
# Regressions
amy.0 ~ SEX_n + AGE
tau.0 ~ a * amy.0 + SEX_n + AGE
WMH   ~ b * amy.0 + SEX_n + AGE
HVR   ~ c * amy.0 + d * tau.0 + e * WMH + SEX_n + AGE
MMSE  ~ f * amy.0 + g * tau.0 + h * WMH + i * HVR + SEX_n + AGE
# Mediation
# Direct effect
dAmy := c
# Indirect effects
iTau  := a * d
iWMH  := b * e
# Total effect
Total := dAmy + iTau + iWMH
# Proportions
pAmy  := dAmy / Total
pTau  := iTau / Total
pWMH  := iWMH / Total
"),
str_glue("
# Regressions
amy.0 ~ SEX_n + AGE
amy.1 ~ SEX_n + AGE
tau.0 ~ SEX_n + AGE
tau.1 ~ SEX_n + AGE
WMH   ~ a0 * amy.0 + a1 * amy.1 + SEX_n + AGE
HVR   ~ b0 * amy.0 + b1 * amy.1 + c0 * tau.0 + c1 * tau.1
  + d * WMH + SEX_n + AGE
MMSE  ~ e0 * amy.0 + e1 * amy.1 + f0 * tau.0 + f1 * tau.1
  + g * WMH + h * HVR + SEX_n + AGE
# Mediation
# Direct effects
dAmy0 := b0
dAmy1 := b1
dTau0 := c0
dTau1 := c1
# Indirect effects
iWMH  := (a0 + a1) * d
# Total effect
Total := dAmy0 + dAmy1 + dTau0 + dTau1 + iWMH
# Proportions
pAmy0 := dAmy0 / Total
pAmy1 := dAmy1 / Total
pTau0 := dTau0 / Total
pTau1 := dTau1 / Total
pWMH  := iWMH / Total
"),
str_glue("
# Regressions
amy.0 ~ SEX_n + AGE
amy.1 ~ SEX_n + AGE
amy.2 ~ SEX_n + AGE
tau.0 ~ SEX_n + AGE
tau.1 ~ a0 * amy.0 + a2 * amy.2 + SEX_n + AGE
tau.2 ~ b * amy.2 + SEX_n + AGE
WMH   ~ c0 * amy.0 + c1 * amy.1 + c2 * amy.2 + SEX_n + AGE
HVR   ~ d0 * amy.0 + d1 * amy.1 + d2 * amy.2
  + e0 * tau.0 + e1 * tau.1 + e2 * tau.2
  + f * WMH + SEX_n + AGE
MMSE  ~ g0 * amy.0 + g1 * amy.1 + g2 * amy.2
  + h0 * tau.0 + h1 * tau.1 + h2 * tau.2
  + i * WMH + j * HVR + SEX_n + AGE
# Mediation
# Direct effect
dAmy0 := d0
dAmy1 := d1
dAmy2 := d2
dTau0 := e0
# Indirect effects
iTau1 := (a0 + a2) * e1
iTau2 := a2 * e2
iWMH  := (c0 + c1 + c2) * f
# Total effect
Total := dAmy0 + dAmy2 + dAmy2 + dTau0 + iTau1 + iTau2 + iWMH
# Proportions
pAmy0 := dAmy0 / Total
pAmy1 := dAmy1 / Total
pAmy2 := dAmy2 / Total
pTau0 := dTau0 / Total
pTau1 := iTau1 / Total
pTau2 := iTau2 / Total
pWMH  := iWMH / Total
"))

hcv.mods    <-
  list(str_glue("
# Regressions
amy.0 ~ SEX_n + AGE
tau.0 ~ a * amy.0 + SEX_n + AGE
WMH   ~ b * amy.0 + SEX_n + AGE
HCv   ~ c * amy.0 + d * tau.0 + e * WMH + SEX_n + AGE
MMSE  ~ f * amy.0 + g * tau.0 + h * WMH + i * HCv + SEX_n + AGE
# Mediation
# Direct effect
dAmy := c
# Indirect effects
iTau  := a * d
iWMH  := b * e
# Total effect
Total := dAmy + iTau + iWMH
# Proportions
pAmy  := dAmy / Total
pTau  := iTau / Total
pWMH  := iWMH / Total
"),
str_glue("
# Regressions
amy.0 ~ SEX_n + AGE
amy.1 ~ SEX_n + AGE
tau.0 ~ SEX_n + AGE
tau.1 ~ SEX_n + AGE
WMH   ~ a0 * amy.0 + a1 * amy.1 + SEX_n + AGE
HCv   ~ b0 * amy.0 + b1 * amy.1 + c0 * tau.0 + c1 * tau.1
  + d * WMH + SEX_n + AGE
MMSE  ~ e0 * amy.0 + e1 * amy.1 + f0 * tau.0 + f1 * tau.1
  + g * WMH + h * HCv + SEX_n + AGE
# Mediation
# Direct effects
dAmy0 := b0
dAmy1 := b1
dTau0 := c0
dTau1 := c1
# Indirect effects
iWMH  := (a0 + a1) * d
# Total effect
Total := dAmy0 + dAmy1 + dTau0 + dTau1 + iWMH
# Proportions
pAmy0 := dAmy0 / Total
pAmy1 := dAmy1 / Total
pTau0 := dTau0 / Total
pTau1 := dTau1 / Total
pWMH  := iWMH / Total
"),
str_glue("
# Regressions
amy.0 ~ SEX_n + AGE
amy.1 ~ SEX_n + AGE
amy.2 ~ SEX_n + AGE
tau.0 ~ SEX_n + AGE
tau.1 ~ a0 * amy.0 + a2 * amy.2 + SEX_n + AGE
tau.2 ~ b * amy.2 + SEX_n + AGE
WMH   ~ c0 * amy.0 + c1 * amy.1 + c2 * amy.2 + SEX_n + AGE
HCv   ~ d0 * amy.0 + d1 * amy.1 + d2 * amy.2
  + e0 * tau.0 + e1 * tau.1 + e2 * tau.2
  + f * WMH + SEX_n + AGE
MMSE  ~ g0 * amy.0 + g1 * amy.1 + g2 * amy.2
  + h0 * tau.0 + h1 * tau.1 + h2 * tau.2
  + i * WMH + j * HCv + SEX_n + AGE
# Mediation
# Direct effect
dAmy0 := d0
dAmy1 := d1
dAmy2 := d2
dTau0 := e0
# Indirect effects
iTau1 := (a0 + a2) * e1
iTau2 := a2 * e2
iWMH  := (c0 + c1 + c2) * f
# Total effect
Total := dAmy0 + dAmy2 + dAmy2 + dTau0 + iTau1 + iTau2 + iWMH
# Proportions
pAmy0 := dAmy0 / Total
pAmy1 := dAmy1 / Total
pAmy2 := dAmy2 / Total
pTau0 := dTau0 / Total
pTau1 := iTau1 / Total
pTau2 := iTau2 / Total
pWMH  := iWMH / Total
"))

### Fit models
## Exploratory models looking for interconnections of Amy/Tau networks
fname       <- here("data/rds/med-mods_nets_raket.rds")

if (!file.exists(fname) | refit_mods) {
  groups    <- DT[, levels(RAKET_group)]
  mod_nets  <- vector("list", length(groups))
  setattr(mod_nets, "names", groups)

  pb <- progress_bar$new(format = "Models | :what [:bar] :current/:total",
                         total = length(groups),
                         clear = FALSE, width = 75)

  for (i in seq_along(groups)) {
    pb$tick(tokens = list(what = sprintf("Networks: %s", groups[i])))
    mod_nets[[i]] <-
      sem(nets.mods[[i]],
          data = DT[RAKET_group == groups[i]],
          estimator = "ML",
          se = "bootstrap",
          bootstrap = 1000)
  }

  write_rds(mod_nets, fname)
} else {
  mod_nets <- read_rds(fname)
}
rm(fname)

## Path plots
labels      <- vector("list", length(mod_nets))
labels[[1]] <- c(AGE = "Age", SEX_n = "Sex", amy.0 = "AB-H", tau.0 = "Tau-H")
labels[[2]] <- c(labels[[1]], amy.1 = "AB-E", tau.1 = "Tau-E")
labels[[3]] <- c(labels[[2]], amy.2 = "AB-L", tau.2 = "Tau-L")

for (i in seq_along(mod_nets)) {
  fname     <- here(sprintf("plots/mediation_paths_nets-%i_%s", i, "full"))
  if (!file.exists(fname) & print_plots) {
    lavaanPlot2(model = mod_nets[[i]], labels = labels[[i]],
                graph_options = list(rankdir = "LR"),
                node_options = list(shape = "box"),
                edge_options = list(color = "grey"),
                coef_labels = T, stand = T, stars = "regress") |>
         embed_plot_pdf(fname)
  }

  fname     <- here(sprintf("plots/mediation_paths_nets-%i_%s", i, "sign"))
  if (!file.exists(fname) & print_plots) {
    coefs   <- extract_coefs(mod_nets[[i]], stand = TRUE) |>
      as.data.table() |>
      {\(x) x[p_val < 0.05]}()
    if (coefs[, .N] == 0) break
    ndf     <- create_nodes(coefs, labels[[i]], list(shape = "box"))
    edf     <- create_edges(coefs, ndf, list(color = "grey"),
                            coef_labels = TRUE, stars = "regress")
    dot     <- convert_graph(ndf, edf, list(rankdir = "LR"))
    lavaanPlot2(gr_viz = dot) |>
      embed_plot_pdf(fname)
    #rm(coefs, ndf, edf, dot)
  }
}

## HVR full models
fname       <- here("data/rds/med-mods_hvr_raket.rds")

if (!file.exists(fname) | refit_mods) {
  groups    <- DT[, levels(RAKET_group)]
  mod_hvr   <- vector("list", length(groups))
  setattr(mod_hvr, "names", groups)

  pb <- progress_bar$new(format = "Models | :what [:bar] :current/:total",
                         total = length(groups),
                         clear = FALSE, width = 75)

  for (i in seq_along(groups)) {
    pb$tick(tokens = list(what = sprintf("HVR: %s", groups[i])))
    mod_hvr[[i]] <-
      sem(hvr.mods[[i]],
          data = DT[RAKET_group == groups[i]],
          estimator = "ML",
          se = "bootstrap",
          bootstrap = 1000)
  }

  write_rds(mod_hvr, fname)
} else {
  mod_hvr   <- read_rds(fname)
}
rm(fname)

## Path plots
for (i in seq_along(mod_hvr)) {
  fname     <- here(sprintf("plots/mediation_paths_hvr-%i_%s", i, "full"))
  if (!file.exists(fname) & print_plots) {
    lavaanPlot2(model = mod_hvr[[i]], labels = labels[[i]],
                graph_options = list(rankdir = "LR"),
                node_options = list(shape = "box"),
                edge_options = list(color = "grey"),
                coef_labels = T, stand = T, stars = "regress") |>
         embed_plot_pdf(fname)
  }

  fname     <- here(sprintf("plots/mediation_paths_hvr-%i_%s", i, "sign"))
  if (!file.exists(fname) & print_plots) {
    coefs   <- extract_coefs(mod_hvr[[i]], stand = TRUE) |>
      as.data.table() |>
      {\(x) x[p_val < 0.05]}()
    if (coefs[, .N] == 0) break
    ndf     <- create_nodes(coefs, labels[[i]], list(shape = "box"))
    edf     <- create_edges(coefs, ndf, list(color = "grey"),
                            coef_labels = TRUE, stars = "regress")
    dot     <- convert_graph(ndf, edf, list(rankdir = "LR"))
    lavaanPlot2(gr_viz = dot) |>
      embed_plot_pdf(fname)
    #rm(coefs, ndf, edf, dot)
  }
}

#fname       <- here("plots/mediation_moca_path.pdf")
#if (!file.exists(fname) | print_plots) {
#labels      <- c(AGE_scan = "Age", SEX_n = "Sex", EDUC = "Education",
                 #HVR_mean_inv = "HC-atrophy", MOCA_score = "MoCA")

  #p_plot    <- lavaanPlot2(model = mod_moca_w.fit, labels = labels,
                            #graph_options = list(rankdir = "LR"),
                            #node_options = list(shape = "box"),
                            #edge_options = list(color = "grey"),
                            #coef_labels = T, stand = T,
                            #stars = "regress")
  #embed_plot_pdf(p_plot, fname)
#}
#rm(fname)

## HCv full models
fname       <- here("data/rds/med-mods_hcv_raket.rds")

if (!file.exists(fname) | refit_mods) {
  groups    <- DT[, levels(RAKET_group)]
  mod_hcv   <- vector("list", length(groups))
  setattr(mod_hcv, "names", groups)

  pb <- progress_bar$new(format = "Models | :what [:bar] :current/:total",
                         total = length(groups),
                         clear = FALSE, width = 75)

  for (i in seq_along(groups)) {
    pb$tick(tokens = list(what = sprintf("HCv: %s", groups[i])))
    mod_hcv[[i]] <-
      sem(hcv.mods[[i]],
          data = DT[RAKET_group == groups[i]],
          estimator = "ML",
          se = "bootstrap",
          bootstrap = 1000)
  }

  write_rds(mod_hcv, fname)
} else {
  mod_hcv   <- read_rds(fname)
}
rm(fname)

## Path plots
for (i in seq_along(mod_hcv)) {
  fname     <- here(sprintf("plots/mediation_paths_hcv-%i_%s", i, "full"))
  if (!file.exists(fname) & print_plots) {
    lavaanPlot2(model = mod_hcv[[i]], labels = labels[[i]],
                graph_options = list(rankdir = "LR"),
                node_options = list(shape = "box"),
                edge_options = list(color = "grey"),
                coef_labels = T, stand = T, stars = "regress") |>
         embed_plot_pdf(fname)
  }

  fname     <- here(sprintf("plots/mediation_paths_hcv-%i_%s", i, "sign"))
  if (!file.exists(fname) & print_plots) {
    coefs   <- extract_coefs(mod_hcv[[i]], stand = TRUE) |>
      as.data.table() |>
      {\(x) x[p_val < 0.05]}()
    if (coefs[, .N] == 0) break
    ndf     <- create_nodes(coefs, labels[[i]], list(shape = "box"))
    edf     <- create_edges(coefs, ndf, list(color = "grey"),
                            coef_labels = TRUE, stars = "regress")
    dot     <- convert_graph(ndf, edf, list(rankdir = "LR"))
    lavaanPlot2(gr_viz = dot) |>
      embed_plot_pdf(fname)
    #rm(coefs, ndf, edf, dot)
  }
}
### Extract Fit measures and paramater estimates
## Names for data cleaning
#msrs_names      <- mod_hvr.fits[[1]] |> fitMeasures() |> names()
#rois_names      <- rois_amy_tau.dt[, paste(LABEL_name, SIDE, sep = "_")]

## HVR model
#mod_hvr.msrs    <- mod_hvr.fits |>
  #lapply(fitMeasures) |>
  #lapply(as.data.table) |>
  #lapply(transpose) |>
  #rbindlist()

#setnames(mod_hvr.msrs, msrs_names)
#mod_hvr.msrs[, ROI := rois_names]

#mod_hvr.est     <- mod_hvr.fits |>
  #lapply(standardizedSolution) |>
  #lapply(as.data.table)

#mod_hvr.est     <- rois_names |>
  #lapply(function (name) mod_hvr.est[[name]][, ROI := name]) |>
  #rbindlist()

#mod_hvr.est     <- mod_hvr.est[op != "~~"]

#mod_moca.msrs   <- mod_moca.fits |>
  #lapply(fitMeasures) |>
  #lapply(as.data.table) |>
  #lapply(transpose) |>
  #rbindlist()

#setnames(mod_moca.msrs, msrs_names)
#mod_hvr.msrs[, ROI := rois_names]

#mod_moca.est    <- mod_moca.fits |>
  #lapply(standardizedSolution) |>
  #lapply(as.data.table)

#mod_moca.est    <- rois_names |>
  #lapply(function (name) mod_moca.est[[name]][, ROI := name]) |>
  #rbindlist()

#mod_moca.est    <- mod_moca.est[op != "~~"]

### Plot standardized estimates
## Use only ROIs in BOTH AMY and TAU important lists


## HVR model
#DT <- mod_hvr.est[lhs %in% c("dAMY", "iTAU") &
                  #ROI %in% rois_amy_tau.dt[LIST == "BOTH",
                                           #paste(LABEL_name, SIDE,
                                                 #sep = "_")]]
#DT[, label := factor(label, levels = c("dAMY", "iTAU"),
                     #labels = c("Direct: Amyloid", "Indirect: Tau"))]

#ordered_rois  <- DT[lhs == "iTAU"][order(est.std), ROI]

#p1 <- DT |>
  #ggplot(aes(x = ROI, y = est.std)) +
  #theme_classic(base_size = 12) +
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  #geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper),
                #width = 0.2, position = position_dodge(0.5)) +
  #geom_point(shape = 21, fill = "white", size = 1.5) +
  #scale_x_discrete(limits = ordered_rois) +
  #labs(x = "CerebrA ROIs", y = "Standardized estimates",
       #title = "Direct & indirect effects on HC atrophy",
       #caption = "Bootstrap CIs: 1000 resamples") +
  #coord_flip() +
  #facet_wrap(vars(label))

##here("plots/mediation_hvr_std-estimates.png") |>
  ##ggsave(p, width = 6, height = 5, units = "in", dpi = 600)

## MoCA model
##mod_moca.est[lhs == "dAMY", label := "Direct: Amyloid"]
##mod_moca.est[lhs == "iTAU", label := "Indirect: Tau"]
##mod_moca.est[lhs == "iHVR", label := "Indirect: HVR"]

##DT <- mod_moca.est[lhs %in% c("dAMY", "iTAU", "iHVR", "Total")]
##DT[, label := factor(label,
                     ##levels = c("Total", "dAMY", "iTAU", "iHVR"),
                     ##labels = c("Total", "Direct: Amyloid",
                                ##"Indirect: Tau", "Indirect: HVR"))]

#DT <- mod_moca.est[lhs %in% c("dAMY", "iTAU") &
                  #ROI %in% rois_amy_tau.dt[LIST == "BOTH",
                                           #paste(LABEL_name, SIDE,
                                                 #sep = "_")]]
#DT[, label := factor(label, levels = c("dAMY", "iTAU"),
                     #labels = c("Direct: Amyloid", "Indirect: Tau"))]

#ordered_rois  <- DT[lhs == "iTAU"][order(-est.std), ROI]

#p2 <- DT |>
  #ggplot(aes(x = ROI, y = est.std)) +
  #theme_classic(base_size = 12) +
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  #geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper),
                #width = 0.2, position = position_dodge(0.5)) +
  #geom_point(shape = 21, fill = "white", size = 1.5) +
  #scale_x_discrete(limits = ordered_rois) +
  #labs(x = "CerebrA ROIs", y = "Standardized estimates",
       #title = "Direct & indirect effects on MoCA scores",
       #caption = "Bootstrap CIs: 1000 resamples") +
  #coord_flip() +
  #facet_wrap(vars(label))

#pp <- p1 + p2
#here("plots/mediation_std-estimates.png") |>
  #ggsave(pp, width = 12, height = 5, units = "in", dpi = 600)


