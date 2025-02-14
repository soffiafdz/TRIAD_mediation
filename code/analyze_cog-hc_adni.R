#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lme4)
#library(effsize)
#library(ggplot2)
#library(GGally)
#library(ggtext)
library(cocor)
library(gt)
library(stargazer)

### INPUT
paths <- list(
  rds = c(
    "adni_age_calculated.rds",
    "adni_hc-hvr_adj.rds",
    "adni_cfa-factors_cog-domains.rds"
  ),
  scripts = c(
    "calculate_age_adni.R",
    "adjust_hc-hvr.R",
    "cfa-cog_adni.R"
  )
) |> Map(
  f = function(Files, Dir) here(Dir, Files),
  c("data/rds", "code")
)

if (file.exists(paths[["rds"]][1])) {
  age.dt <- readRDS(paths[["rds"]][1])
} else {
  sprintf("Sourcing: %s", paths[["scripts"]][1])
  source(paths[["scripts"]][1])
}

if (file.exists(paths[["rds"]][2])) {
  hc_hvr.dt <- readRDS(paths[["rds"]][2])
} else {
  sprintf("Sourcing: %s", paths[["scripts"]][2])
  source(paths[["scripts"]][2])
}

if (file.exists(paths[["rds"]][3])) {
  cog_lat.dt <- readRDS(paths[["rds"]][3])
} else {
  sprintf("Sourcing: %s", paths[["scripts"]][3])
  source(paths[["scripts"]][3])
}

rm(paths)

### Data CLEANING
## Age & Sex
library(ADNIMERGE)
data(adnimerge)
setDT(adnimerge)
age_sex.dt    <- adnimerge[, .(SEX = PTGENDER), PTID] |>
  unique() |>
  merge(age.dt) |>
  setkey(PTID, VISCODE)
rm(adnimerge, age.dt)

## Labels
adjs          <- c("Unadjusted", "Proportions", "Residuals")
dxs           <- c("CH", "MCI", "AD")
rois          <- c(
  ICC = "Intracranial volume",
  HC  = "Hippocampus",
  VC  = "Ventricles",
  HVR = "HC-to-VC ratio"
)
cog_factors <- c(
  "Executive function", "Language", "Memory", "Visospatial function"
)
names(cog_factors) <- names(cog_lat.dt[, -c("PTID", "VISCODE")])

## Average sides for HC/HVR & Make long format
hc_hvr.dt <- hc_hvr.dt[
    ,
    # Average sides for all ROIs, except ICC
    lapply(.SD, mean),
    by = .(PTID, VISCODE, DX, ICC, ADJ),
    .SDcols = names(rois)[-1]
  ][
    , DX := factor(DX, levels = c("CN", "Dementia", "MCI"), labels = dxs)
  ][
    # Rename ADJ for plotting
    , ADJ := factor(ADJ, labels = adjs)
  ][
    # ICC is always unadjusted, so delete duplicate values
    !adjs[1], on = "ADJ", ICC := NA
  ] |>
  # Rename ROIs for plotting
  setnames(names(rois), rois) |>
  # Aggregate ROIs
  melt(measure = rois, variable = "ROI", value = "VAL", na.rm = TRUE) |>
  # Key resulting DT
  setkey(PTID, VISCODE)

## Make lat_cog long format
cog_lat.dt  <- cog_lat.dt |>
  melt(id = 1:2, variable = "COGDOM", value = "SCORE") |>
  na.omit() |>
  setkey(PTID, VISCODE)

### Data PROCESSING
cors.dt <- list()
cor.cols <- list(
  R = "estimate.cor",
  CI_l = "conf.int1",
  CI_h = "conf.int2",
  STAT = "statistic.t",
  P_val = "p.value"
)

mods.lst <- list()
for (fct in names(cog_factors)) {
  mods.lst[[fct]] <- list()
  for (roi in names(rois)) {
    mods.lst[[fct]][[roi]] <- list()
    for (adj in adjs) {
      DT <- hc_hvr.dt[
        adj, on = "ADJ"
      ][
        rois[roi], on = "ROI",
        .(DX, VAL),
        .(PTID, VISCODE)
      ][
        cog_lat.dt[fct, on = "COGDOM"],
        nomatch = NULL
      ]
      if (DT[, .N] > 0) {
        # Correlations
        for (dx in dxs) {
          subDT <- DT[dx, on = "DX"]
          # N of subsample
          n_sub <- subDT[!duplicated(PTID), .N, ]
          n_vis <- subDT[, .N, ]
          if (n_vis > 0) {
            cor.res <- subDT[, cor.test(VAL, SCORE)] |> unlist()
            sublist <- list(
              DX = dx,
              COGDOM = fct,
              ADJ = adj,
              ROI = rois[roi],
              N_subs = n_sub,
              N_vis = n_vis
            )
            for (i in seq_along(cor.cols)) {
              sublist[[names(cor.cols[i])]] <- cor.res[[cor.cols[[i]]]]
            }
            setDT(sublist)
            cors.dt <- rbind(cors.dt, sublist)
            rm(i, cor.res, sublist)
          }
        }
        ## MLM models
        # Std coefficients based on earliest timepoint of Ctrl group
        mods.lst[[fct]][[roi]][[adj]] <- DT |>
          setkey(PTID, VISCODE) |>
          merge(age_sex.dt) |>
          {function(DT) {
            ctrl_stats <- DT[
              "CH",
              on = "DX",
              .SD[which.min(EXAMDATE)],
              PTID
            ][
              , .(
                VAL_mean = mean(VAL),
                VAL_sd = sd(VAL),
                SCORE_mean = mean(SCORE),
                SCORE_sd = sd(SCORE)
              )
            ]
            DT[
              , let(
                VAL = (VAL - ctrl_stats$VAL_mean) / ctrl_stats$VAL_sd,
                SCORE = (SCORE - ctrl_stats$SCORE_mean) / ctrl_stats$SCORE_sd
              )
            ]
          }} () |>
          lmer(formula = SCORE ~ AGE + SEX + DX + VAL + (1|PTID)) |>
          suppressWarnings()
      }
    }
  }
}

cors.dt[
  , (names(cor.cols)) := lapply(.SD, as.numeric),
  .SDcols = names(cor.cols)
][
  , P_adj := p.adjust(P_val, method = "bonferroni")
]
rm(dx, fct, roi, adj, DT, subDT, n_vis, n_sub, cor.cols)

## Test correlation differences
# Need correlation between ROIs
# N is different by CogDomain, so it needs to be done by group
intercors.dt <- hc_hvr.dt[
  (ROI %like% "Hip|Vent" & ADJ == "Residuals") |
  (ROI %like% "ratio" & ADJ == "Proportions"),
  -"ADJ"
][
  cog_lat.dt[, -"SCORE"], # Remove Scores since they are not needed
  nomatch = NULL
] |>
  dcast(... ~ ROI, value.var = "VAL") |>
  {
    function(DT)
    DT[
      ,
      {
        cormatrix <- cor(.SD)
        list(
          HC_VC = cormatrix[rois[2], rois[3]],
          HC_HVR = cormatrix[rois[2], rois[4]],
          VC_HVR = cormatrix[rois[3], rois[4]]
        )
      },
      keyby = .(COGDOM, DX),
      .SDcols = rois[-1] # There is no ICC exploration
    ]
  } ()

# Merge
cordiff.dt <- cors.dt[
  (ROI %like% "Hip|Vent" & ADJ == "Residuals") |
  (ROI %like% "ratio" & ADJ == "Proportions"),
  .(
    N = N_vis,
    ROI = factor(ROI, levels = rois[-1], labels = names(rois[-1])),
    R
  ),
  keyby = .(COGDOM, DX)
] |>
  dcast(... ~ ROI, value.var = "R") |>
  merge(intercors.dt)

DT <- list()
for (i in seq_len(cordiff.dt[, .N])) {
  for (XY in c("HC_VC", "HC_HVR", "VC_HVR")) {
    X <- strsplit(XY, split = "_")[[1]][1]
    Y <- strsplit(XY, split = "_")[[1]][2]
    cordiff <- cocor.dep.groups.overlap(
      r.jk = cordiff.dt[i, get(X)],
      r.jh = cordiff.dt[i, get(Y)],
      r.kh = cordiff.dt[i, get(XY)],
      n    = cordiff.dt[i, N]
    ) |> get.cocor.results()
    DT <- rbind(
      DT,
      cordiff.dt[
        i,
        .(
          COGDOM,
          DX,
          N,
          COMPARISON = XY,
          Z = cordiff$hittner2003$statistic,
          p = cordiff$hittner2003$p.value,
          CI_l = cordiff$zou2007$conf.int[1],
          CI_h = cordiff$zou2007$conf.int[2]
        )
      ]
    )
  }
}
cordiff.dt <- copy(DT)
cordiff.dt[, p_adj := p.adjust(p, method = "bonferroni")]
rm(i, XY, X, Y, cordiff, DT)

### Print tables
# TODO: Condition this by a starting-script-constant
for (fct in names(cog_factors)) {
  DT <- cors.dt[
    !(ADJ == adjs[3] & ROI == rois["HVR"]) ## Remove HVR (Residuals)
  #][
    #P_adj < 0.05 ## Keep only significant correlations
  ][
    fct, on = "COGDOM"
  ]

  # Remove adjustment labels for HVR & ICC
  DT[rois[c(1,4)], on = "ROI", ADJ := NA]

  for (all_adjs in c(TRUE, FALSE)) {
    fpath <- sprintf(
      "tables/cors_cog_%s%s.html",
      tolower(fct),
      ifelse(all_adjs, "_all_adjs", "")
    ) |> here()

    if (all_adjs) {
      subDT <- copy(DT)
    } else {
      subDT <- DT[
        !ADJ %in% c("Unadjusted", "Proportions")
      ][
        !rois["ICC"],
        on = "ROI"
      ]
      subDT[, ADJ := NA]
      browser()
    }

    subDT[order(-abs(R)), -c("COGDOM", "DX", "N_subs", "N_vis", "P_val")] |>
    gt() |>
    tab_header(
      title = "Correlation of cognition with volumetry",
      subtitle = cog_factors[fct]
    ) |>
    fmt_number(decimals = 3) |>
    cols_merge(columns = c("ROI", "ADJ"), pattern = "{1}<< ({2})>>") |>
    cols_merge(columns = c("R", "CI_l", "CI_h"), pattern = "{1} ({2}, {3})") |>
    cols_label(
      ROI = md(ifelse(all_adjs, "**ROI (adj.)**", "**ROI**")),
      R = md("***r***"),
      STAT = md("**T**"),
      P_adj = md("***p* adj.**")
    ) |>
    tab_row_group(
      label = sprintf(
        "Dementia (N: %s, timepoints: %s)",
        subDT["AD", on = "DX", N_subs] |> unique() |> format(big.mark = ","),
        subDT["AD", on = "DX", N_vis] |> unique() |> format(big.mark = ",")
      ),
      rows = subDT[order(-abs(R)), which(DX == "AD")]
    ) |>
    tab_row_group(
      label = sprintf(
        "Mild cognitive impairment (N: %s, timepoints: %s)",
        subDT["MCI", on = "DX", N_subs] |> unique() |> format(big.mark = ","),
        subDT["MCI", on = "DX", N_vis] |> unique() |> format(big.mark = ",")
      ),
      rows = subDT[order(-abs(R)), which(DX == "MCI")]
    ) |>
    tab_row_group(
      label = sprintf(
        "Cognitively healthy (N: %s, timepoints: %s)",
        subDT["CH", on = "DX", N_subs] |> unique() |> format(big.mark = ","),
        subDT["CH", on = "DX", N_vis] |> unique() |> format(big.mark = ",")
      ),
      rows = subDT[order(-abs(R)), which(DX == "CH")]
    ) |>
    gtsave(fpath)
  }

  ### Mixed-effect models
  ## Compare adjustment methods
  # ICC + HC
  fpath <- fct |>
    tolower() |>
    sprintf(fmt = "tables/mlm_%s_adjs_icc-hc.html") |>
    here()

  mods.lst[[fct]][c("ICC", "HC")] |>
    stargazer(
      out = fpath,
      intercept.bottom = FALSE,
      type = "html",
      #single.row = TRUE,
      title = "MLM: Cognition ~ Head-size & Hippocampus",
      column.labels = rois[c("ICC", "HC")],
      column.separate = c(1, 3),
      dep.var.labels.include = FALSE,
      dep.var.caption = cog_factors[fct],
      covariate.labels = c(
        "Intercept", "Age (y)", "Sex (M)", "Dx (MCI)", "Dx (AD)", "Volume (CC)"
      ),
      notes = paste(
        "Head-size adjustments methods were:",
        "(2) Unadjusted, (3) Proportions, (4) Residuals"
      )
    )

  # ICC + VC
  fpath <- fct |>
    tolower() |>
    sprintf(fmt = "tables/mlm_%s_adjs_icc-vc.html") |>
    here()

  mods.lst[[fct]][c("ICC", "VC")] |>
    stargazer(
      out = fpath,
      intercept.bottom = FALSE,
      type = "html",
      #single.row = TRUE,
      title = "MLM: Cognition ~ Head-size & Ventricles",
      column.labels = rois[c("ICC", "HC")],
      column.separate = c(1, 3),
      dep.var.labels.include = FALSE,
      dep.var.caption = cog_factors[fct],
      covariate.labels = c(
        "Intercept", "Age (y)", "Sex (M)", "Dx (MCI)", "Dx (AD)", "Volume (CC)"
      ),
      notes = paste(
        "Head-size adjustments methods were:",
        "(2) Unadjusted, (3) Proportions, (4) Residuals"
      )
    )

  fpath <- fct |>
    tolower() |>
    sprintf(fmt = "tables/mlm_%s_adjs_icc-hvr.html") |>
    here()

  mods.lst[[fct]]["HVR"] |>
    stargazer(
      out = fpath,
      intercept.bottom = FALSE,
      type = "html",
      #single.row = TRUE,
      title = "MLM: Cognition ~ Hippocampus-to-Ventricle ratio",
      column.labels = c("Regular", "Robust"),
      model.numbers = FALSE,
      dep.var.labels.include = FALSE,
      dep.var.caption = cog_factors[fct],
      covariate.labels = c(
        "Intercept", "Age (y)", "Sex (M)", "Dx (MCI)", "Dx (AD)", "HVR"
      ),
      notes = "Robust HVR was calculated using residuals-adjusted volumes."
    )

  ## ICC, {HC, VC}-residuals & HVR
  fpath <- fct |>
    tolower() |>
    sprintf(fmt = "tables/mlm_%s_icc-hc-vc-hvr.html") |>
    here()
  # There is a bug in stargazer, models need to be saved in an object
  submods <- c(
    mods.lst[[fct]][["ICC"]],
    mods.lst[[fct]][["HC"]][["Residuals"]],
    mods.lst[[fct]][["VC"]][["Residuals"]],
    mods.lst[[fct]][["HVR"]][["Proportions"]]
  )

  stargazer(
    submods,
    out = fpath,
    intercept.bottom = FALSE,
    type = "html",
    #single.row = TRUE,
    title = "MLM: Cognition ~ Medio-temporal volumetry",
    column.labels = rois,
    model.numbers = FALSE,
    dep.var.labels.include = FALSE,
    dep.var.caption = cog_factors[fct],
    covariate.labels = c(
      "Intercept", "Age (y)", "Sex (M)", "Dx (MCI)", "Dx (AD)", "Volume/HVR"
    ),
    notes = paste(
      "Hippocampal and Ventricular volumes were adjusted for head-size",
      "using the residuals method"
    )
  )
}
