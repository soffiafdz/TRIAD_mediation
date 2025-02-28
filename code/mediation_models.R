#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lavaan)
library(lavaanPlot)

## Print plots
print_plots <- TRUE

# Load baseline data
fpath <- here("data/rds/triad.rds")
if (file.exists(fpath)) {
  triad.dt <- readRDS(fpath)
} else {
  here("code/demographics.R") |> source()
}
rm(fpath)

## Mediation analysis: Imaging
# Keep only subjects with full imaging
# Remove youth and AD
triad.dt <- triad.dt[
  !is.na(AMYLOID)
][
  !is.na(TAU_braak1)
][
  !DX_clean %in% c("Young", "Other", "AD"),
  .(
    ## Subj ID
    PTID,
    ### Group
    #DX = factor(DX, levels = c("CN", "MCI")),
    ## Exogenous variable
    AMYLOID,
    ## Endogenous variables
    TAU = rowSums(.SD), # Sum Braak Stages
    # 1 - HVR (average for both sides)
    HVR_inv = 1 - (HVR_l + HVR_r) / 2,
    MOCA = MOCA_score,
    ## Covariates
    AGE = AGE_scan,
    EDUC,
    # Convert Sex to dummy variable
    SEX = as.numeric(SEX) - 1,  ## M: 1 & F: 0
    # APOE status (N of alleles)
    APOE_n,
    ## Sex & APOE4 moderation
    SEX_APOE = APOE_n * (as.numeric(SEX) - 1)
  ),
  .SDcols = TAU_braak1:TAU_braak6
]

## Labels:
labels.v <- c(
  AGE = "Age",
  SEX = "Sex",
  AMY = "Amyloid (PET)",
  TAU = "Tau (PET)",
  APOE_n = "APOE4",
  EDUC = "Education",
  SEX_APOE = "Sex*APOE4",
  HVR_inv = "HC-atrophy",
  MOCA = "MoCA"
)


### Original AAIC model
## AMYLOID -> TAU (sum) -> 1-HVR_mean -> COG
mod_orig.lst <- list()
mod_orig.lst[["MODEL"]] <- '
  # Regression equations for serial moderation
  AMYLOID ~ c1*AGE + c2*SEX + c3*APOE_n

  TAU ~ a*AMYLOID + t1*AGE + t2*SEX + t3*APOE_n

  HVR_inv ~ b*AMYLOID + c*TAU + h1*AGE + h2*SEX + h3*APOE_n

  MOCA ~ d*AMYLOID + e*TAU + f*HVR_inv + m1*AGE + m2*EDUC + m3*SEX

  # Total effect: direct, simple indirect, and serial effects
  Total := d + (a*e) + (b*f) + (a*c*f)

  # Direct effect of Amyloid on MoCA
  deAMY := d
  propAMY := deAMY / Total

  # Indirect effect via TAU only
  ieTAU := a * e
  propTAU := ieTAU / Total

  # Indirect effect via HVR only
  ieHVR := b * f
  propHVR := ieHVR / Total

  # Serial mediation Amyloid -> Tau -> HVR -> MoCA
  ieSerial := a * c * f
  propSerial := ieSerial / Total

  # Combined indirect effect
  ieTotal := (a*e) + (b*f) + (a*c*f)
  propTotal := ieTotal/Total
'

fname <- here('data/rds/mediation_aaic2025.rds')
if (file.exists(fname)) {
  mod_orig.lst[["FIT"]] <- readRDS(fname)
} else {
  mod_orig.lst[["FIT"]] <- sem(
    mod_orig.lst[["MODEL"]],
    data = triad.dt,
    cluster = "PTID",
    estimator = "ML",
    se = "robust.cluster",
    bootstrap = 2000
  )
  saveRDS(mod_orig.lst[["FIT"]], fname)
}

rm(fname)

if (print_plots) {
  coefs <- extract_coefs(mod_orig.lst[["FIT"]], stand = TRUE) |> setDT()
  fpaths <- c("skel", "reg") |>
    sprintf(fmt = "plots/med_aaic2025_%s.pdf") |>
    here()

  ## Structural model (skeleton):
  coefs["~", on = "op"] |>
  {\(coefs) {
    ndf <- create_nodes(coefs, labels.v, NULL)
    edf <- create_edges(coefs, ndf, list(color = "grey"), coef_labels = FALSE)
    dot <- convert_graph(ndf, edf, list(rankdir = "LR"))
    lavaanPlot2(gr_viz = dot)
    }}() |> embed_plot_pdf(fpaths[1])

  ## Structural model (Significant paths):
  coefs["~", on = "op"]["" != stars] |>
  {\(coefs) {
    ndf <- create_nodes(coefs, labels.v, NULL)
    edf <- create_edges(
      coefs,
      ndf,
      list(color = "grey"),
      coef_labels = TRUE,
      stars = "regress"
    )
    dot <- convert_graph(ndf, edf, list(rankdir = "LR"))
    lavaanPlot2(gr_viz = dot)
    }}() |> embed_plot_pdf(fpaths[2])
}

### Updated model with moderation
## AMYLOID -> TAU (sum) -> 1-HVR_mean -> COG
mod.lst <- list()
mod.lst[["MODEL"]] <- '
  # Regression equations for serial moderation
  AMYLOID ~ c1*AGE + c2*SEX + c3*APOE_n + c4*SEX_APOE

  TAU ~ a*AMYLOID +
    t1*AGE + t2*SEX + t3*APOE_n + t4*SEX_APOE

  HVR_inv ~ b*AMYLOID + c*TAU +
    h1*AGE + h2*SEX + h3*APOE_n + h4*SEX_APOE

  MOCA ~ d*AMYLOID + e*TAU + f*HVR_inv +
    m1*AGE + m2*EDUC + m3*SEX

  # Total effect: direct, simple indirect, and serial effects
  Total := d + (a*e) + (b*f) + (a*c*f)

  # Direct effect of Amyloid on MoCA
  deAMY := d
  propAMY := deAMY / Total

  # Indirect effect via TAU only
  ieTAU := a * e
  propTAU := ieTAU / Total

  # Indirect effect via HVR only
  ieHVR := b * f
  propHVR := ieHVR / Total

  # Serial mediation Amyloid -> Tau -> HVR -> MoCA
  ieSerial := a * c * f
  propSerial := ieSerial / Total

  # Combined indirect effect
  ieTotal := (a*e) + (b*f) + (a*c*f)
  propTotal := ieTotal/Total
'

fname <- here('data/rds/mediation_aaic2025_moderation.rds')
if (file.exists(fname)) {
  mod.lst[["FIT"]] <- readRDS(fname)
} else {
  mod.lst[["FIT"]] <- sem(
    mod.lst[["MODEL"]],
    data = triad.dt,
    cluster = "PTID",
    estimator = "ML",
    se = "robust.cluster",
    bootstrap = 2000
  )
  saveRDS(mod.lst[["FIT"]], fname)
}

rm(fname)

if (print_plots) {
  coefs <- extract_coefs(mod.lst[["FIT"]], stand = TRUE) |> setDT()
  fpaths <- c("skel", "reg") |>
    sprintf(fmt = "plots/med_aaic2025_moderation_%s.pdf") |>
    here()

  ## Structural model (skeleton):
  coefs["~", on = "op"] |>
  {\(coefs) {
    ndf <- create_nodes(coefs, labels.v, NULL)
    edf <- create_edges(coefs, ndf, list(color = "grey"), coef_labels = FALSE)
    dot <- convert_graph(ndf, edf, list(rankdir = "LR"))
    lavaanPlot2(gr_viz = dot)
    }}() |> embed_plot_pdf(fpaths[1])

  ## Structural model (Significant paths):
  coefs["~", on = "op"]["" != stars] |>
  {\(coefs) {
    ndf <- create_nodes(coefs, labels.v, NULL)
    edf <- create_edges(
      coefs,
      ndf,
      list(color = "grey"),
      coef_labels = TRUE,
      stars = "regress"
    )
    dot <- convert_graph(ndf, edf, list(rankdir = "LR"))
    lavaanPlot2(gr_viz = dot)
    }}() |> embed_plot_pdf(fpaths[2])
}
