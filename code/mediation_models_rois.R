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

# Wide → Long format
DT <- melt(DT, measure = patterns("amy", "tau"),
           variable = "NET", value = c("AMY", "TAU"))

# Keep only Disease with their respective networks
DT <- rbindlist(list(DT[RAKET_group %like% "Healthy" & NET == 1],
                     DT[RAKET_group %like% "Early" & NET == 2],
                     DT[RAKET_group %like% "Late" & NET == 3]))

## Model definitions
hvr.mod     <- str_glue("
# Regressions
AMY ~ SEX_n + AGE
TAU ~ a * AMY + SEX_n + AGE
WMH   ~ b * AMY + SEX_n + AGE
HVR   ~ c * AMY + d * TAU + e * WMH + SEX_n + AGE
MMSE  ~ f * AMY + g * TAU + h * WMH + i * HVR + SEX_n + AGE
# Mediation
# Direct effect
dAmy := f
# Indirect effects
iTau  := a * g
iWMH  := b * h
iHVR  := (c + d + e) * i
# Total effect
Total := dAmy + iTau + iWMH + iHVR
# Proportions
pAmy  := dAmy / Total
pTau  := iTau / Total
pWMH  := iWMH / Total
pHVR  := iHVR / Total
")

hcv.mod     <-
  str_glue("
# Regressions
AMY ~ SEX_n + AGE
TAU ~ a * AMY + SEX_n + AGE
WMH   ~ b * AMY + SEX_n + AGE
HCv   ~ c * AMY + d * TAU + e * WMH + SEX_n + AGE
MMSE  ~ f * AMY + g * TAU + h * WMH + i * HCv + SEX_n + AGE
# Mediation
# Direct effect
dAmy := f
# Indirect effects
iTau  := a * g
iWMH  := b * h
iHCv  := (c + d + e) * i
# Total effect
Total := dAmy + iTau + iWMH + iHCv
# Proportions
pAmy  := dAmy / Total
pTau  := iTau / Total
pWMH  := iWMH / Total
pHCv  := iHCv / Total
")

### Fit models
## Path plots
## HVR full models
fname       <- here("data/rds/med-mod_hvr_raket_ug.rds")

if (!file.exists(fname) | refit_mods) {
  mod_hvr <-
    sem(hvr.mod,
        data = DT,
        estimator = "ML",
        se = "bootstrap",
        bootstrap = 1000)

  write_rds(mod_hvr, fname)
} else {
  mod_hvr   <- read_rds(fname)
}
rm(fname)

## Path plots
fname     <- here("plots/mediation_paths_hvr_ug_full.pdf")
if (!file.exists(fname) & print_plots) {
  lavaanPlot2(model = mod_hvr, labels = c(SEX_n = "SEX"),
              graph_options = list(rankdir = "LR"),
              node_options = list(shape = "box"),
              edge_options = list(color = "grey"),
              coef_labels = T, stand = T, stars = "regress") |>
       embed_plot_pdf(fname)
}

fname     <- here("plots/mediation_paths_hvr_ug_eskeleton.pdf")
if (!file.exists(fname) & print_plots) {
  lavaanPlot2(model = mod_hvr, labels = c(SEX_n = "SEX"),
              graph_options = list(rankdir = "LR"),
              node_options = list(shape = "box"),
              edge_options = list(color = "grey"),
              coef_labels = F) |>
       embed_plot_pdf(fname)
}

fname     <- here("plots/mediation_paths_hvr_ug_sign.pdf")
if (!file.exists(fname) & print_plots) {
  coefs   <- extract_coefs(mod_hvr, stand = TRUE) |>
    as.data.table() |>
    {\(x) x[p_val < 0.05]}()
  if (coefs[, .N] == 0) next
  ndf     <- create_nodes(coefs, c(SEX_n = "SEX"), list(shape = "box"))
  edf     <- create_edges(coefs, ndf, list(color = "grey"),
                          coef_labels = TRUE, stars = "regress")
  dot     <- convert_graph(ndf, edf, list(rankdir = "LR"))
  lavaanPlot2(gr_viz = dot) |>
    embed_plot_pdf(fname)
  #rm(coefs, ndf, edf, dot)
}

## HCv full models
fname       <- here("data/rds/med-mod_hcv_raket_ug.rds")
if (!file.exists(fname) | refit_mods) {
  mod_hcv <-
    sem(hcv.mod,
        data = DT,
        estimator = "ML",
        se = "bootstrap",
        bootstrap = 1000)

  write_rds(mod_hcv, fname)
} else {
  mod_hcv   <- read_rds(fname)
}
rm(fname)

## Path plots
fname     <- here("plots/mediation_paths_hcv_ug_full.pdf")
if (!file.exists(fname) & print_plots) {
  lavaanPlot2(model = mod_hcv, labels = c(SEX_n = "SEX"),
              graph_options = list(rankdir = "LR"),
              node_options = list(shape = "box"),
              edge_options = list(color = "grey"),
              coef_labels = T, stand = T, stars = "regress") |>
       embed_plot_pdf(fname)
}

fname     <- here("plots/mediation_paths_hcv_ug_sign.pdf")
if (!file.exists(fname) & print_plots) {
  coefs   <- extract_coefs(mod_hcv, stand = TRUE) |>
    as.data.table() |>
    {\(x) x[p_val < 0.05]}()
  if (coefs[, .N] == 0) next
  ndf     <- create_nodes(coefs, c(SEX_n = "SEX"), list(shape = "box"))
  edf     <- create_edges(coefs, ndf, list(color = "grey"),
                          coef_labels = TRUE, stars = "regress")
  dot     <- convert_graph(ndf, edf, list(rankdir = "LR"))
  lavaanPlot2(gr_viz = dot) |>
    embed_plot_pdf(fname)
  #rm(coefs, ndf, edf, dot)
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


