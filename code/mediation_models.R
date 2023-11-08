#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lavaan)
library(lavaanPlot)

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

## Amyloid -> TAU -> HC (HVR)
simple1.mod <- '
  # Latent variables
  HVR =~ HVR_l + HVR_r
  TAU =~ TAU_braak1 + TAU_braak2 + TAU_braak3 +
    TAU_braak4 + TAU_braak5 + TAU_braak6
  # Direct effect
  HVR ~ c * AMYLOID + SEX_n + AGE_scan + APOE_n
  # Mediator
  TAU ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  HVR ~ b * TAU
  # Indirect effect
  indirect := a * b
  # Total effect
  total := c + (a*b)
'

simple1.fit <- sem(simple1.mod, data = triad.dt, estimator = "ML")

## amyloid -> TAU -> HC (HCvol)
simple2.mod <- '
  # Latent variables
  HC =~ HCv_l + HCv_r
  TAU =~ TAU_braak1 + TAU_braak2 + TAU_braak3 +
    TAU_braak4 + TAU_braak5 + TAU_braak6
  # Direct effect
  HC ~ c * AMYLOID + SEX_n + AGE_scan + APOE_n
  # Mediator
  TAU ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  HC ~ b * TAU
  # Indirect effect
  indirect := a * b
  # Total effect
  total := c + (a*b)
'

simple2.fit <- sem(simple2.mod, data = triad.dt, estimator = "ML")

## Neuropsych
# Keep only subjects with full neuropsy evaluation (dependent variable)
# MoCA scores and three RAVLT subscores
#triad.dt <- triad.dt[!is.na(RAVLT_intro) & !is.na(MOCA_score)]
triad.dt <- triad.dt[!is.na(MOCA_score)]
triad.dt[, MOCA_scaled := scale(MOCA_score)]

## Amyloid -> HC (HVR) -> MoCA
simple3.mod <- '
  # Latent variables
  HVR =~ HVR_l + HVR_r
  # Direct effect
  MOCA_score ~ c * AMYLOID + SEX_n + AGE_scan + APOE_n + EDUC
  # Mediator
  HVR ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  MOCA_score ~ b * HVR
  # Indirect effect
  indirect := a * b
  # Total effect
  total := c + (a*b)
'

simple3.fit <- sem(simple3.mod, data = triad.dt, estimator = "ML")

## Include memory: Regressions are not significant
triad.dt <- triad.dt[!is.na(RAVLT_intro)]
triad.dt[, `:=`(RAVLT_intro = scale(RAVLT_intro),
                   RAVLT_raw = scale(RAVLT_raw),
                   RAVLT_rep = scale(RAVLT_rep))]

## Amyloid -> HC (HVR) -> Cognition
simple4.mod <- '
  # Latent variables
  HVR =~ HVR_l + HVR_r
  MEM =~ RAVLT_rep
  # Direct effect
  MEM ~ c * AMYLOID + SEX_n + AGE_scan + APOE_n + EDUC
  # Mediator
  HVR ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  MEM ~ b * HVR
  # Indirect effect
  indirect := a * b
  # Total effect
  total := c + (a*b)
'

simple4.fit <- sem(simple4.mod, data = triad.dt, estimator = "ML")

## Amyloid -> HC (HVR) -> Cognition
simple5.mod <- '
  # Latent variables
  HVR =~ HVR_l + HVR_r
  MEM =~ RAVLT_rep + RAVLT_raw + RAVLT_intro
  NEURO =~ MEM + MOCA_score
  # Direct effect
  NEURO ~ c * AMYLOID + SEX_n + AGE_scan + APOE_n + EDUC
  # Mediator
  HVR ~ a * AMYLOID + SEX_n + AGE_scan + APOE_n
  NEURO ~ b * HVR
  # Indirect effect
  indirect := a * b
  # Total effect
  total := c + (a*b)
'

simple5.fit <- sem(simple5.mod, data = triad.dt, estimator = "ML")

# Serial Mediation
# AMYLOID -> TAU -> HVR -> COG
