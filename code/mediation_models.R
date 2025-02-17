#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lavaan)
library(lavaanPlot)

## Print plots
print_plots <- TRUE

# Load baseline data
fpath       <- here("data/rds/triad.rds")
if (file.exists(fpath)) {
  triad.dt  <- readRDS(fpath)
} else {
  here("code/demographics.R") |> source()
}

## Mediation analysis: Imaging
# Keep only subjects with full imaging
# Remove youth and AD
triad.dt    <- triad.dt[
  !is.na(AMYLOID)
][
  !is.na(TAU_braak1)
][
  !DX_clean %in% c("Young", "Other", "AD")
]

## Cleaning
triad.dt[
  ,
  let(
    # Convert Sex to dummy variable
    SEX_n = as.numeric(SEX) - 1,
    DX = factor(DX, levels = c("CN", "MCI")),
    # 1 - HVR (average for both sides)
    HVR_lr = 1 - HVR_l,
    HVR_rr = 1 - HVR_r,
    HVR_mean_inv = 1 - (HVR_l + HVR_r) / 2,
    # Sum Braak Stages
    TAU = rowSums(.SD)
  ),
  .SDcols = TAU_braak1:TAU_braak6
]

## Labels:
labels.lst  <- list(
  COV = c(AGE_scan = "Age", SEX_n = "Sex", APOE_n = "APOE4"),
  HCV = c(HCv_l = "Left", HCv_r = "Right"),
  HVR = c(
    HVR = "HC-atrophy",
    HVR_lr = "1-HVR (Left)",
    HVR_rr = "1-HVR (Right)"
  ),
  TAU = c(
    TAU_braak1 = "Braak1",
    TAU_braak2 = "Braak2",
    TAU_braak3 = "Braak3",
    TAU_braak4 = "Braak4",
    TAU_braak5 = "Braak5",
    TAU_braak6 = "Braak6"
  ),
  HVR2 = c(HVR_mean_inv = "HC-atrophy"),
  MOCA = c(MOCA_score = "MoCA"),
  MEM = c(
    RAVLT_rep = "RAVLT (rep)",
    RAVLT_intro = "RAVLT (intro)",
    RAVLT_raw = "RAVLT (raw)")
)


## AMYLOID -> TAU (sum) -> 1-HVR_mean -> COG
serial7.mod <- '
  # Regressions
  AMYLOID ~ SEX_n + AGE_scan + APOE_n
  TAU ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan + APOE_n
  MOCA_score ~ d * AMYLOID + e * TAU + f * HVR_mean_inv + SEX_n + AGE_scan + APOE_n + EDUC
  # Direct effect
  deAMY := d
  # Indirect effects
  ieTAU := a * e
  ieHVR := b * f
  SerialMed := a * c * f
  # Total effect
  Total := deAMY + ieTAU + ieHVR + SerialMed
  # Proportions
  propAMY := deAMY / Total
  propTAU := ieTAU / Total
  propHVR := ieHVR / Total
  propSerial := SerialMed / Total
'

fname <- here('data/rds/mediation_aaic2025.rds')
if (file.exists(fname)) {
  serial7.fit <- readRDS(fname)
} else {
  serial7.fit <- sem(
    serial7.mod,
    data = triad_cog.dt,
    estimator = "ML",
    se = "bootstrap",
    bootstrap = 10000
  )
  saveRDS(serial7.fit, fname)
}

rm(fname)

if (print_plots) {
  ## TODO: Change to labels.lst
  labels      <- c(
    labels.lst[["COV"]],
    labels.lst[["HVR2"]],
    labels.lst[["MOCA"]]
  )

  p_ser7      <- lavaanPlot2(
    model = serial7.fit,
    labels = labels,
    graph_options = list(rankdir = "LR"),
    node_options = list(shape = "box"),
    edge_options = list(color = "grey"),
    coef_labels = T,
    stand = T,
    stars = "regress"
  )

  here("data/derivatives/med_aaic2025.pdf") |> embed_plot_pdf(plot = p_ser7)
}

### AMYLOID -> TAU (sum) -> 1-HVR_mean -> COG
#serial7.mod <- '
  ## Regressions
  #AMYLOID ~ SEX_n + AGE_scan + APOE_n
  #TAU ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  #HVR_mean_inv ~ b * AMYLOID + c * TAU + SEX_n + AGE_scan + APOE_n
  #MOCA_score ~ d * AMYLOID + e * TAU + f * HVR_mean_inv + SEX_n + AGE_scan + APOE_n + EDUC
  ## Direct effect
  #deAMY := d
  #propAMY := deAMY / Total
  ## Indirect effects
  #ieTAU := a * e
  #propTAU := ieTAU / Total
  #ieHVR := b * f
  #propHVR := ieHVR / Total
  #SerialMed := a * c * f
  #propSerial := SerialMed / Total
  ## Total effect
  #Total := deAMY + ieTAU + ieHVR + SerialMed
#'

#fname <- here('data/rds/mediation_serial7.rds')
#if (file.exists(fname)) {
  #serial7.fit <- readRDS(fname)
#} else {
  #serial7.fit <- sem(
    #serial7.mod,
    #data = triad_cog.dt,
    #estimator = "ML",
    #se = "bootstrap",
    #bootstrap = 10000
  #)
  #saveRDS(serial7.fit, fname)
#}
#rm(fname)
#if (print_plots) {
  ### TODO: Change to labels.lst
  #labels      <- c(
    #labels.lst[["COV"]],
    #labels.lst[["HVR2"]],
    #labels.lst[["MOCA"]]
  #)

  #p_ser7      <- lavaanPlot2(
    #model = serial7.fit,
    #labels = labels,
    #graph_options = list(rankdir = "LR"),
    #node_options = list(shape = "box"),
    #edge_options = list(color = "grey"),
    #coef_labels = T,
    #stand = T,
    #stars = "regress"
  #)

  #here("data/derivatives/serial7.pdf") |> embed_plot_pdf(plot = p_ser7)
#}
