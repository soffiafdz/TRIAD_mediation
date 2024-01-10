#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
library(lavaan)
library(lavaanPlot)

## Print plots ### Needs to be done outside renv
print_plots <- FALSE
#print_plots <- TRUE

# Load baseline data
fpath       <- here("data/rds/triad.rds")
if (file.exists(fpath)) {
  triad.dt  <- readr::read_rds(fpath)
} else {
  here("code/demographics.R") |> source()
}

# Convert Sex to dummy variable
triad.dt[, SEX_n := as.numeric(SEX) - 1]

## Mediation analysis: Imaging
# Keep only subjects with full imaging
triad.dt    <- triad.dt[!is.na(AMYLOID) & !is.na(TAU_braak1)]

# Remove youth and AD
triad.dt    <-triad.dt[!DX_clean %in% c("Young", "Other", "AD")]
triad.dt[, DX := factor(DX, levels = c("CN", "MCI"))]

# 1 - HVR (average for both sides)
triad.dt[, `:=`(HVR_lr = 1 - HVR_l, HVR_rr = 1 - HVR_r)]
triad.dt[, `:=`(HVR_mean_inv = 1 - (HVR_l + HVR_r) / 2)]

# Sum Braak Stages
triad.dt[, TAU := TAU_braak1 + TAU_braak2 + TAU_braak3 +
         TAU_braak4 + TAU_braak5 + TAU_braak6]

## Labels:
labels_cov  <- c(AGE_scan = "Age", SEX_n = "Sex", APOE_n = "APOE4")
labels_hcv  <- c(HCv_l = "Left", HCv_r = "Right")
labels_hvr  <- c(HVR = "HC-atrophy",
                 HVR_lr = "1-HVR (Left)",
                 HVR_rr = "1-HVR (Right)")
labels_tau  <- c(TAU_braak1 = "Braak1",
                 TAU_braak2 = "Braak2",
                 TAU_braak3 = "Braak3",
                 TAU_braak4 = "Braak4",
                 TAU_braak5 = "Braak5",
                 TAU_braak6 = "Braak6")
labels_hvr2 <- c(HVR_mean_inv = "HC-atrophy")
labels_moca <- c(MOCA_score = "MoCA")
labels_mem  <- c(RAVLT_rep = "RAVLT (rep)",
                 RAVLT_intro = "RAVLT (intro)",
                 RAVLT_raw = "RAVLT (raw)")

## Amyloid -> TAU -> HC (HVR)
simple1.mod <- '
  # Latent variables
  HVR =~ HVR_l + HVR_rr
  TAU =~ TAU_braak1 + TAU_braak2 + TAU_braak3 +
    TAU_braak4 + TAU_braak5 + TAU_braak6
  # Regressions
  AMYLOID ~ SEX_n + AGE_scan + APOE_n
  TAU ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  HVR ~ b * TAU + c * AMYLOID + SEX_n + AGE_scan + APOE_n
  # Indirect effect
  indirect := a * b
  # Total effect
  total := c + (a*b)
'
#simple1.fit <- sem(simple1.mod, data = triad.dt, estimator = "ML")
#if (print_plots) {
  #labels      <- c(labels_cov, labels_hvr, labels_tau)
  #p_simp1     <- lavaanPlot2(model = simple1.fit, labels = labels,
                             #graph_options = list(rankdir = "LR"),
                             #node_options = list(shape = "box"),
                             #edge_options = list(color = "grey"),
                             #coef_labels = T, stand = T,
                             ##stars = c("regress", "latent")
                             #)
  #embed_plot_pdf(p_simp1, here("data/derivatives/simple1.pdf"))
#}

## amyloid -> TAU -> HC (HCvol)
simple2.mod <- '
  # Latent variables
  HCvol =~ HCv_l + HCv_r
  TAU =~ TAU_braak1 + TAU_braak2 + TAU_braak3 +
    TAU_braak4 + TAU_braak5 + TAU_braak6
  # Regressions
  AMYLOID ~ SEX_n + AGE_scan + APOE_n
  TAU ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  HCvol ~ b * TAU + c * AMYLOID + SEX_n + AGE_scan + APOE_n
  # Indirect effect
  indirect := a * b
  # Total effect
  total := c + (a*b)
'

#simple2.fit <- sem(simple2.mod, data = triad.dt, estimator = "ML")
#if (print_plots) {
  #labels      <- c(labels_cov, labels_hcv, labels_tau)
  #p_simp2     <- lavaanPlot2(model = simple2.fit, labels = labels,
                             #graph_options = list(rankdir = "LR"),
                             #node_options = list(shape = "box"),
                             #edge_options = list(color = "grey"),
                             #coef_labels = T, stand = T,
                             #stars = c("regress", "latent"))
  #embed_plot_pdf(p_simp2, here("data/derivatives/simple2.pdf"))
#}

## Neuropsych
# Keep only subjects with full neuropsy evaluation (dependent variable)
# MoCA scores and three RAVLT subscores
#triad.dt <- triad.dt[!is.na(RAVLT_intro) & !is.na(MOCA_score)]
triad_cog.dt <- triad.dt[!is.na(MOCA_score)]
triad_cog.dt[, MOCA_scaled := scale(MOCA_score)]

## Amyloid -> HC (HVR) -> MoCA
simple3.mod <- '
  # Latent variables
  HVR =~ HVR_lr + HVR_rr
  # Regressions
  AMYLOID ~ SEX_n + AGE_scan + APOE_n
  HVR ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  MOCA_score ~ b * HVR + c * AMYLOID + SEX_n + AGE_scan + APOE_n + EDUC
  # Indirect effect
  indirect := a * b
  # Total effect
  total := c + (a*b)
'

#simple3.fit <- sem(simple3.mod, data = triad_cog.dt, estimator = "ML")
#if (print_plots) {
  #labels      <- c(labels_cov, labels_hvr, labels_moca)
  #p_simp3     <- lavaanPlot2(model = simple3.fit, labels = labels,
                            #graph_options = list(rankdir = "LR"),
                            #node_options = list(shape = "box"),
                            #edge_options = list(color = "grey"),
                            #coef_labels = T, stand = T,
                            #stars = c("regress", "latent"))
  #embed_plot_pdf(p_simp3, here("data/derivatives/simple3.pdf"))
#}

## Include memory: Regressions are not significant
triad_mem.dt <- triad_cog.dt[!is.na(RAVLT_intro)]
triad_mem.dt[, `:=`(RAVLT_intro = scale(RAVLT_intro),
                    RAVLT_raw = scale(RAVLT_raw),
                    RAVLT_rep = scale(RAVLT_rep))]

## Amyloid -> HC (HVR) -> Memory
simple4.mod <- '
  # Latent variables
  HVR =~ HVR_lr + HVR_rr
  MEMORY =~ RAVLT_rep + RAVLT_raw + RAVLT_intro
  # Regressions
  AMYLOID ~ SEX_n + AGE_scan + APOE_n
  HVR ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  MEMORY ~ b * HVR + c * AMYLOID + SEX_n + AGE_scan + APOE_n + EDUC
  # Indirect effect
  indirect := a * b
  # Total effect
  total := c + (a*b)
'

#simple4.fit <- sem(simple4.mod, data = triad_mem.dt, estimator = "ML")
#if (print_plots) {
  #labels      <- c(labels_cov, labels_hvr, labels_mem)
  #p_simp4     <- lavaanPlot2(model = simple4.fit, labels = labels,
                            #graph_options = list(rankdir = "LR"),
                            #node_options = list(shape = "box"),
                            #edge_options = list(color = "grey"),
                            #coef_labels = T, stand = T,
                            #stars = c("regress", "latent"))
  #embed_plot_pdf(p_simp4, here("data/derivatives/simple4.pdf"))
#}

# AMY -> TAU1 -> TAU2 -> TAU3 -> TAU4 -> TAU5 -> TAU6
simple5.mod <- '
  # Latent variables
  TAU =~ TAU_braak1 + TAU_braak2 + TAU_braak3 +
    TAU_braak4 + TAU_braak5 + TAU_braak6
  # Regressions
  AMYLOID ~ SEX_n + AGE_scan + APOE_n
  TAU_braak1 ~ a1 * AMYLOID + SEX_n + AGE_scan + APOE_n
  TAU_braak2 ~ a2 * AMYLOID + b1 * TAU_braak1 + SEX_n + AGE_scan + APOE_n
  TAU_braak3 ~ a3 * AMYLOID + b2 * TAU_braak1 + c1 * TAU_braak2 +
    SEX_n + AGE_scan + APOE_n
  TAU_braak4 ~ a4 * AMYLOID + b3 * TAU_braak1 + c2 * TAU_braak2 +
    d1 * TAU_braak3 + SEX_n + AGE_scan + APOE_n
  TAU_braak5 ~ a5 * AMYLOID + b4 * TAU_braak1 + c3 * TAU_braak2 +
    d2 * TAU_braak3 + e1 * TAU_braak4 + SEX_n + AGE_scan + APOE_n
  TAU_braak6 ~ a6 * AMYLOID + b5 * TAU_braak1 + c4 * TAU_braak2 +
    d3 * TAU_braak3 + e2 * TAU_braak4 + f1 *TAU_braak5 +
    SEX_n + AGE_scan + APOE_n
  # Indirect effect
  #indirect := a * b
  # Total effect
  #total := c + (a*b)
'
#simple5.fit <- sem(simple5.mod, data = triad.dt, estimator = "ML")
#if (print_plots) {
  #labels      <- c(labels_cov, labels_tau)
  #p_simp5     <- lavaanPlot2(model = simple5.fit, labels = labels,
                             #graph_options = list(rankdir = "LR"),
                             #node_options = list(shape = "box"),
                             #edge_options = list(color = "grey"),
                             #coef_labels = T, stand = T,
                             #stars = c("regress", "latent")
                             #)
  #embed_plot_pdf(p_simp1, here("data/derivatives/simple5.pdf"))
#}

# Serial Mediation
# AMYLOID -> TAU -> HVR -> COG
serial1.mod <- '
  # Latent variables
  HVR =~ HVR_lr + HVR_rr
  TAU =~ TAU_braak1 + TAU_braak2 + TAU_braak3 +
    TAU_braak4 + TAU_braak5 + TAU_braak6
  # Regressions
  AMYLOID ~ SEX_n + AGE_scan + APOE_n
  TAU ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  HVR ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan + APOE_n
  MOCA_score ~ d * AMYLOID + e * TAU + f * HVR + SEX_n + AGE_scan + APOE_n + EDUC
  # Total effect
  Total := d + (a * e) + (b * f)
  # Indirect effects
  ieAMY := d
  propAMY := ieAMY / Total
  ieTAU := a * e
  propTAU := ieTAU / Total
  ieHVR := b * f
  propHVR := ieHVR / Total
  ieTotal := (a * e) + (b * f)
  propTotal := ieTotal/ Total

'

fname <- here('data/rds/mediation_serial1.rds')
if (file.exists(fname)) {
  serial1.fit <- read_rds(fname)
} else {
  serial1.fit <- sem(serial1.mod, data = triad_cog.dt, estimator = "ML")
  #serial1.fit <- sem(serial1.mod, data = triad_cog.dt, estimator = "ML",
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(serial1.fit, fname)
}
rm(fname)
#serial1g.fit <- sem(serial1.mod, data = triad_cog.dt, estimator = "ML",
                    #group = "DX")
if (print_plots) {
  labels      <- c(labels_cov, labels_hvr, labels_tau, labels_moca)
  p_ser1      <- lavaanPlot2(model = serial1.fit, labels = labels,
                            graph_options = list(rankdir = "LR"),
                            node_options = list(shape = "box"),
                            edge_options = list(color = "grey"),
                            coef_labels = T, stand = T,
                            stars = c("regress", "latent"))
  embed_plot_pdf(p_ser1, here("data/derivatives/serial1.pdf"))
}

# AMYLOID -> TAU -> HCv -> COG
serial2.mod <- '
  # Latent variables
  HCvol =~ HCv_l + HCv_r
  TAU =~ TAU_braak1 + TAU_braak2 + TAU_braak3 +
    TAU_braak4 + TAU_braak5 + TAU_braak6
  # Regressions
  AMYLOID ~ SEX_n + AGE_scan + APOE_n
  TAU ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  HCvol ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan + APOE_n
  MOCA_score ~ d * AMYLOID + e * TAU + f * HCvol
    + SEX_n + AGE_scan + APOE_n + EDUC
  # Indirect effects
  ieTAU := a * e
  ieHCvol := b * f
  ieTotal := (a * e) + (b * f)
  # Total effect
  Total := d + (a * e) + (b * f)
'

#serial2.fit <- sem(serial2.mod, data = triad_cog.dt, estimator = "ML")
#if (print_plots) {
  #labels      <- c(labels_cov, labels_hcv, labels_tau, labels_moca)
  #p_ser2      <- lavaanPlot2(model = serial2.fit, labels = labels,
                            #graph_options = list(rankdir = "LR"),
                            #node_options = list(shape = "box"),
                            #edge_options = list(color = "grey"),
                            #coef_labels = T, stand = T,
                            #stars = c("regress", "latent"))
  #embed_plot_pdf(p_ser2, here("data/derivatives/serial2.pdf"))
#}

# AMYLOID -> TAU -> HCv -> MEM + COG
serial3.mod <- '
  # Latent variables
  HVR =~ HVR_lr + HVR_rr
  TAU =~ TAU_braak1 + TAU_braak2 + TAU_braak3 +
    TAU_braak4 + TAU_braak5 + TAU_braak6
  MEMORY =~ RAVLT_rep + RAVLT_raw + RAVLT_intro
  # Regressions
  AMYLOID ~ SEX_n + AGE_scan + APOE_n
  TAU ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  HVR ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan + APOE_n
  #MOCA_score ~ d * AMYLOID + e * TAU + f * HVR + SEX_n + AGE_scan + APOE_n + EDUC
  MEMORY ~ g * AMYLOID + h * TAU + i * HVR + SEX_n + AGE_scan + APOE_n + EDUC
  ### MOCA
  ## Indirect effects
  #ieTAU_cog := a * e
  #ieHVR_cog := b * f
  #ieTotal_cog := (a * e) + (b * f)
  ## Total effect
  #Total_cog := d + (a * e) + (b * f)
  ## MEMORY
  ieTAU_mem := a * h
  ieHVR_mem := b * i
  ieTotal_mem := (a * h) + (b * i)
  # Total effect
  Total_mem := g + (a * h) + (b * i)
  # Covariances
  #MEMORY ~~ MOCA_score
'

#serial3.fit <- sem(serial3.mod, data = triad_mem.dt, estimator = "ML")
#if (print_plots) {
  #labels      <- c(labels_cov, labels_hvr, labels_tau, labels_mem)
  #p_ser3      <- lavaanPlot2(model = serial3.fit, labels = labels,
                            #graph_options = list(rankdir = "LR"),
                            #node_options = list(shape = "box"),
                            #edge_options = list(color = "grey"),
                            #coef_labels = T, stand = T,
                            #stars = c("regress", "latent"))
  #embed_plot_pdf(p_ser3, here("data/derivatives/serial3.pdf"))
#}

serial4.mod <- '
  # Latent variables
  HVR =~ HVR_lr + HVR_rr
  TAU =~ TAU_braak1 + TAU_braak2 + TAU_braak3 +
    TAU_braak4 + TAU_braak5 + TAU_braak6
  # Regressions
  AMYLOID ~ SEX_n + AGE_scan + APOE_n
  TAU ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  HVR ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan + APOE_n
  MMSE ~ d * AMYLOID + e * TAU + f * HVR + SEX_n + AGE_scan + APOE_n + EDUC
  # Total effect
  Total := d + (a * e) + (b * f)
  # Indirect effects
  ieAMY := d
  propAMY := ieAMY / Total
  ieTAU := a * e
  propTAU := ieTAU / Total
  ieHVR := b * f
  propHVR := ieHVR / Total
  ieTotal := (a * e) + (b * f)
  propTotal := ieTotal/ Total
'

fname <- here('data/rds/mediation_serial4.rds')
if (file.exists(fname)) {
  serial4.fit <- read_rds(fname)
} else {
  serial4.fit <- sem(serial4.mod, data = triad_cog.dt, estimator = "ML")
                     #se = "bootstrap", bootstrap = 5000)
  write_rds(serial4.fit, fname)
}
rm(fname)

#if (print_plots) {
  #labels      <- c(labels_cov, labels_hvr, labels_tau, labels_moca)
  #p_ser4      <- lavaanPlot2(model = serial4.fit, labels = labels,
                            #graph_options = list(rankdir = "LR"),
                            #node_options = list(shape = "box"),
                            #edge_options = list(color = "grey"),
                            #coef_labels = T, stand = T,
                            #stars = c("regress", "latent"))
  #embed_plot_pdf(p_ser4, here("data/derivatives/serial4.pdf"))
#}

serial5.mod <- '
  # Latent variables
  HVR =~ HVR_lr + HVR_rr
  TAU =~ TAU_braak1 + TAU_braak2 + TAU_braak3 +
    TAU_braak4 + TAU_braak5 + TAU_braak6
  COG =~ MOCA_score + MMSE
  # Regressions
  AMYLOID ~ SEX_n + AGE_scan + APOE_n
  TAU ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  HVR ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan + APOE_n
  COG ~ d * AMYLOID + e * TAU + f * HVR + SEX_n + AGE_scan + APOE_n + EDUC
  # Total effect
  Total := d + (a * e) + (b * f)
  # Indirect effects
  ieAMY := d
  propAMY := ieAMY / Total
  ieTAU := a * e
  propTAU := ieTAU / Total
  ieHVR := b * f
  propHVR := ieHVR / Total
  ieTotal := (a * e) + (b * f)
  propTotal := ieTotal/ Total

'

fname <- here('data/rds/mediation_serial5.rds')
if (file.exists(fname)) {
  serial5.fit <- read_rds(fname)
} else {
  serial5.fit <- sem(serial5.mod, data = triad_cog.dt, estimator = "ML")
                     #se = "bootstrap", bootstrap = 5000)
  write_rds(serial5.fit, fname)
}
rm(fname)

#if (print_plots) {
  #labels      <- c(labels_cov, labels_hvr, labels_tau, labels_moca)
  #p_ser5      <- lavaanPlot2(model = serial5.fit, labels = labels,
                            #graph_options = list(rankdir = "LR"),
                            #node_options = list(shape = "box"),
                            #edge_options = list(color = "grey"),
                            #coef_labels = T, stand = T,
                            #stars = c("regress", "latent"))
  #embed_plot_pdf(p_ser5, here("data/derivatives/serial5.pdf"))
#}

# AMYLOID -> TAU -> WMH -> HVR -> COG
serial6.mod <- '
  # Latent variables
  HVR =~ HVR_lr + HVR_rr
  TAU =~ TAU_braak1 + TAU_braak2 + TAU_braak3 +
    TAU_braak4 + TAU_braak5 + TAU_braak6
  # Regressions
  AMYLOID ~ SEX_n + AGE_scan + APOE_n
  TAU ~ a1 * AMYLOID + SEX_n + AGE_scan + APOE_n
  WMH ~ a2 * AMYLOID + b1 * TAU + SEX_n + AGE_scan + APOE_n
  HVR ~ a3 * AMYLOID + b2 * TAU + c1 * WMH + SEX_n + AGE_scan + APOE_n
  MOCA_score ~ a * AMYLOID + b * TAU + c * WMH + d * HVR +
    SEX_n + AGE_scan + APOE_n + EDUC
  # Direct effect
  # Indirect effects
  ieAMY := a
  ieTAU := a1 * b
  ieWMH := a2 * c
  ieHVR := a3 * d
  ieTotal := ieTAU + ieWMH + ieHVR
  # Total effect
  Total := ieAMY + ieTotal
  # Proportion mediations
  propAMY := ieAMY / Total
  propTAU := ieTAU / Total
  protWMH := ieWMH / Total
  propHVR := ieHVR / Total
  propTotal := ieTotal/ Total
'

fname <- here('data/rds/mediation_serial6.rds')
if (file.exists(fname)) {
  serial6.fit <- read_rds(fname)
} else {
  serial6.fit <- sem(serial6.mod, data = triad_cog.dt, estimator = "ML")
  #serial1.fit <- sem(serial1.mod, data = triad_cog.dt, estimator = "ML",
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(serial6.fit, fname)
}
rm(fname)
#serial1g.fit <- sem(serial1.mod, data = triad_cog.dt, estimator = "ML",
                    #group = "DX")
if (print_plots) {
  labels      <- c(labels_cov, labels_hvr, labels_tau, labels_moca)
  p_ser6      <- lavaanPlot2(model = serial6.fit, labels = labels,
                            graph_options = list(rankdir = "LR"),
                            node_options = list(shape = "box"),
                            edge_options = list(color = "grey"),
                            coef_labels = T, stand = T,
                            stars = c("regress", "latent"))
  embed_plot_pdf(p_ser6, here("data/derivatives/serial6.pdf"))
}

# AMYLOID -> TAU (sum) -> 1-HVR_mean -> COG
serial7.mod <- '
  # Regressions
  AMYLOID ~ SEX_n + AGE_scan + APOE_n
  TAU ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan + APOE_n
  MOCA_score ~ d * AMYLOID + e * TAU + f * HVR_mean_inv + SEX_n + AGE_scan + APOE_n + EDUC
  # Total effect
  Total := d + (a * e) + (b * f)
  # Indirect effects
  ieAMY := d
  propAMY := ieAMY / Total
  ieTAU := a * e
  propTAU := ieTAU / Total
  ieHVR := b * f
  propHVR := ieHVR / Total
  ieTotal := (a * e) + (b * f)
  propTotal := ieTotal/ Total

'

fname <- here('data/rds/mediation_serial7.rds')
if (file.exists(fname)) {
  serial7.fit <- read_rds(fname)
} else {
  serial7.fit <- sem(serial7.mod, data = triad_cog.dt, estimator = "ML")
  #serial1.fit <- sem(serial1.mod, data = triad_cog.dt, estimator = "ML",
                     #se = "bootstrap", bootstrap = 10000)
  write_rds(serial7.fit, fname)
}
rm(fname)
#serial1g.fit <- sem(serial1.mod, data = triad_cog.dt, estimator = "ML",
                    #group = "DX")
if (print_plots) {
  labels      <- c(labels_cov, labels_hvr2, labels_moca)
  p_ser7      <- lavaanPlot2(model = serial7.fit, labels = labels,
                            graph_options = list(rankdir = "LR"),
                            node_options = list(shape = "box"),
                            edge_options = list(color = "grey"),
                            coef_labels = T, stand = T,
                            stars = "regress")
  embed_plot_pdf(p_ser7, here("data/derivatives/serial7.pdf"))
}

