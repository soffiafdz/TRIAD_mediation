#!/usr/bin/env Rscript

library(here)
library(data.table)
library(ggplot2)
library(ggtext)
library(ggsignif)
# library(ggridges)
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
age_sex.dt <- adnimerge[, .(SEX = PTGENDER), PTID] |>
  unique() |>
  merge(age.dt) |>
  setkey(PTID, VISCODE)
rm(adnimerge, age.dt)

## Labels
adjs <- c("Unadjusted", "Proportions", "Residuals")
dxs <- c("CH", "MCI", "AD")
rois <- c(
  ICC = "Intracranial volume",
  HC  = "Hippocampus",
  VC  = "Ventricles",
  HVR = "HC-to-Vent. ratio"
)
cog_factors <- c(
  "Executive function", "Language", "Memory", "Visospatial function"
)
names(cog_factors) <- names(cog_lat.dt[, -c("PTID", "VISCODE", "DATE")])

## Average sides for HC/HVR & Make long format
hc_hvr.dt <- hc_hvr.dt[
  ,
  # Average sides for all ROIs, except ICC
  lapply(.SD, mean),
  by = .(PTID, VISCODE, DX, ICC, ADJ),
  .SDcols = names(rois)[-1]
][
  , DX := factor(DX, levels = c("CN", "MCI", "Dementia"), labels = dxs)
][
  # Rename ADJ for plotting
  , ADJ := factor(ADJ, labels = adjs)
][
  # ICC is always unadjusted, so delete duplicate values
  !adjs[1],
  on = "ADJ", ICC := NA
] |>
  # Rename ROIs for plotting
  setnames(names(rois), rois) |>
  # Aggregate ROIs
  melt(measure = rois, variable = "ROI", value = "VAL", na.rm = TRUE) |>
  # Key resulting DT
  setkey(PTID, VISCODE)

DT <- cog_lat.dt[
  age_sex.dt,
  nomatch = NULL
][
  hc_hvr.dt,
  on = .(PTID, VISCODE), nomatch = NULL
] |>
  melt(
    measure = c("EXECFUN", "LANGUAGE", "MEMORY", "VISOSPATIALFUN"),
    variable = "COGDOM",
    value = "SCORE"
  ) |>
  na.omit() |>
  setkey(PTID, VISCODE)

## Demographics
setkey(DT, COGDOM)
for (cog in names(cog_factors)) {
  col_order <- c("DX", "SEX", "AGE", "SESSN", cog, names(rois)[-1])
  subDT <- DT[
    cog
  ][
    (ROI %like% "Hip|Ventric" & ADJ == "Residuals") |
      (ROI %like% "ratio" & ADJ == "Proportions")
  ]
  subDT[, ROI := factor(ROI, labels = names(rois)[-1])] |> invisible()
  subDT <- dcast(subDT[, -"ADJ"], ... ~ ROI, value.var = "VAL")
  blDT <- subDT[, .SD[which.min(AGE)], PTID] |> setkey(DX)
  # Standardize according to Controls
  blDT[, (cog) := (SCORE - blDT["CH", mean(SCORE)]) / blDT["CH", sd(SCORE)]] |>
    invisible()
  # Number of sessions
  blDT[subDT[, .(SESSN = .N), PTID], on = "PTID", ..col_order] |>
    setnames(cog, paste(cog_factors[cog], "(Z)")) |>
    tbl_summary(
      by = DX,
      label = list(
        SESSN ~ "Visits (N)",
        AGE ~ "Age (y)",
        SEX ~ "Sex (F)",
        HC ~ sprintf("%s (CC)", rois["HC"]),
        VC ~ sprintf("Lat. %s (CC)", rois["VC"]),
        HVR ~ sprintf("%s (0-1)", rois["HVR"])
      ),
      statistic = list(
        all_continuous() ~ "{mean} ({sd})",
        all_continuous2() ~ c("{max}", "{sum}")
      ),
      value = list(SEX = "Female"),
      type = list(SESSN ~ "continuous2", SEX ~ "dichotomous")
    ) |>
    modify_header(label ~ "**Variable**") |>
    modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Clinical label**") |>
    as_gt() |>
    gtsave(sprintf("demog_%s.html", tolower(cog)), here("tables"))
}


### Data PROCESSING
cors.dt <- cor.cols <- list()
cor.cols$R <- "estimate.cor"
cor.cols$CI_l <- "conf.int1"
cor.cols$CI_h <- "conf.int2"
cor.cols$STAT <- "statistic.t"
cor.cols$P_val <- "p.value"

setkey(DT, DX, COGDOM, ROI, ADJ)
sort.dt <- CJ(DX = dxs, COGDOM = names(cog_factors), ROI = rois, ADJ = adjs)
for (i in sort.dt[, seq(.N)]) {
  subDT <- DT[sort.dt[i]] |> na.omit()
  if (subDT[, .N == 0]) next
  # Correlations
  cor.res <- subDT[, cor.test(VAL, SCORE)] |> unlist()
  sublist <- list(
    DX = sort.dt[i, DX],
    COGDOM = sort.dt[i, COGDOM],
    ADJ = sort.dt[i, ADJ],
    ROI = sort.dt[i, ROI],
    N_subs = subDT[!duplicated(PTID), .N],
    N_vis = subDT[, .N]
  )
  for (j in seq_along(cor.cols)) {
    sublist[[names(cor.cols[j])]] <- as.numeric(cor.res[[cor.cols[[j]]]])
  }
  setDT(sublist)
  cors.dt <- rbind(cors.dt, sublist)
}
rm(i, j, cor.cols, cor.res, subDT, sublist)

# Multiple comparisons correction
cors.dt[, P_adj := p.adjust(P_val, method = "bonferroni")] |> invisible()

## Test correlation differences
# Need correlation between ROIs
# N is different by CogDomain, so it needs to be done by group
intercors.dt <- DT[
  (ROI %like% "Hip|Ventric" & ADJ == "Residuals") |
    (ROI %like% "ratio" & ADJ == "Proportions"),
  .(PTID, VISCODE, DX, COGDOM, ROI, VAL)
] |>
  dcast(... ~ ROI, value.var = "VAL") |>
  {
    \(subDT)
    subDT[
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
  }()

# Merge
cordiff.dt <- cors.dt[
  (ROI %like% "Hip|Ventr" & ADJ == "Residuals") |
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

## Use absolutes for correlation strength
cols <- names(cordiff.dt)[-1:-3]
cordiff.dt[, (cols) := lapply(.SD, abs), .SDcols = cols] |> invisible()
cortest.dt <- data.table()
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
    cortest.dt <- rbind(
      cortest.dt,
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
cortest.dt[, p_adj := p.adjust(p, method = "bonferroni")] |> invisible()
rm(i, XY, X, Y, cordiff)

### Print tables
# TODO: Condition this by a starting-script-constant
setkey(cors.dt, COGDOM)
for (fct in names(cog_factors)) {
  # Remove HVR residuals
  subDT <- cors.dt[fct][!(ADJ == adjs[3] & ROI == rois["HVR"])]
  # Remove adjustment labels for HVR & ICC
  subDT[rois[c(1, 4)], on = "ROI", ADJ := NA]

  for (all_adjs in c(TRUE, FALSE)) {
    fpath <- sprintf(
      "tables/cors_cog_%s%s.html",
      tolower(fct),
      ifelse(all_adjs, "_all_adjs", "")
    ) |> here()

    if (all_adjs) {
      subset.dt <- copy(subDT)
    } else {
      subset.dt <- subDT[
        !ADJ %in% c("Unadjusted", "Proportions")
      ][
        !rois["ICC"],
        on = "ROI"
      ]
      subset.dt[, ADJ := NA]
    }

    subset.dt[order(-abs(R)), -c("COGDOM", "DX", "N_subs", "N_vis", "P_val")] |>
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
          subset.dt[DX == "AD", N_subs] |> unique() |> format(big.mark = ","),
          subset.dt[DX == "AD", N_vis] |> unique() |> format(big.mark = ",")
        ),
        rows = subset.dt[order(-abs(R)), which(DX == "AD")]
      ) |>
      tab_row_group(
        label = sprintf(
          "Mild cognitive impairment (N: %s, timepoints: %s)",
          subset.dt[DX == "MCI", N_subs] |> unique() |> format(big.mark = ","),
          subset.dt[DX == "MCI", N_vis] |> unique() |> format(big.mark = ",")
        ),
        rows = subset.dt[order(-abs(R)), which(DX == "MCI")]
      ) |>
      tab_row_group(
        label = sprintf(
          "Cognitively healthy (N: %s, timepoints: %s)",
          subset.dt[DX == "CH", N_subs] |> unique() |> format(big.mark = ","),
          subset.dt[DX == "CH", N_vis] |> unique() |> format(big.mark = ",")
        ),
        rows = subset.dt[order(-abs(R)), which(DX == "CH")]
      ) |>
      gtsave(fpath)
  }
}

### Plot
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

subDT <- cors.dt[
  (ROI %like% "Hip|Ventr" & ADJ == "Residuals") |
    (ROI %like% "ratio" & ADJ == "Proportions")
][
  !"VISOSPATIALFUN",
  on = "COGDOM", .(DX, COGDOM, ROI, R, CI_l, CI_h, P_adj)
]

subDT[, ROI := factor(ROI, levels = rev(rois), labels = names(rev(rois)))] |>
  invisible()

subDT[, DX := factor(
  DX,
  levels = c("CH", "MCI", "AD") # ,
  # labels = c("Unimpaired", "MCI", "Alzheimer's")
)] |> invisible()

subDT[, COGDOM := factor(
  COGDOM,
  labels = c("Executive Function", "Language", "Memory")
)] |> invisible()

subDT[, LABEL := fcase(
  P_adj < 0.001, "***",
  P_adj < 0.01, "**",
  P_adj < 0.05, "*"
)] |> invisible()

subDT["VC", on = "ROI", let(R = R * -1, CI_h = CI_l * -1, CI_l = CI_h * -1)] |>
  invisible()

subDT2 <- cortest.dt[!COGDOM %like% "V" & p_adj < 0.05, c(1:2, 4, 9)]
subDT2[, DX := factor(DX, levels = c("CH", "MCI", "AD"))] |> invisible()
subDT2[, COGDOM := factor(
  COGDOM,
  labels = c("Executive Function", "Language", "Memory")
)] |> invisible()
subDT2[, c("XMIN", "XMAX") := tstrsplit(COMPARISON, "_")] |> invisible()
subDT2[, LABEL := fcase(
  p_adj < 0.001, "***",
  p_adj < 0.01, "**",
  p_adj < 0.05, "*"
)] |> invisible()
subDT2 <- subDT[, .(Y = max(CI_h)), .(DX, COGDOM)][subDT2, on = .(DX, COGDOM)]
subDT2[order(-COMPARISON), id := 1:.N, .(DX, COGDOM)] |> invisible()
subDT2[, Y := Y + .01 + (id - 1)^1.75 * .015] |> invisible()
subDT2[id == 2 & XMAX == "HVR", Y := Y + .02] |> invisible()
subDT2[, c("id", "p_adj", "COMPARISON") := NULL]

p <- ggplot(subDT, aes(ROI, R)) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(size = 11),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    # plot.caption = element_text(size = 8),
    legend.position = "bottom"
  ) +
  facet_grid(rows = vars(DX), cols = vars(COGDOM), scales = "free_x") +
  geom_errorbar(
    aes(ymin = CI_l, ymax = CI_h, colour = ROI),
    width = 0.2
  ) +
  geom_point(
    aes(colour = ROI),
    shape = 21,
    fill = "white",
    size = 2
    # stroke = 0.3
  ) +
  geom_text(
    aes(label = LABEL, y = R, colour = ROI),
    size = 2.5,
    vjust = -.8
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    alpha = .5,
    colour = cbPalette[1]
  ) +
  geom_signif(
    aes(
      xmin = XMIN,
      xmax = XMAX,
      annotations = LABEL,
      y_position = Y
    ),
    subDT2,
    manual = TRUE,
    colour = cbPalette[1],
    textsize = 3,
    tip_length = 0,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = cbPalette[c(2:3, 8)],
    labels = c("HC-to-Ventricle Ratio", "Lateral Ventricles", "Hippocampus"),
    guide = guide_legend(reverse = TRUE)
  ) +
  scale_y_continuous(expand = expansion(mult = 0.075)) +
  labs(y = "Pearson's R") +
  coord_flip()

fpath <- here("plots/aaic2025_cors.png")
ggsave(fpath, p, width = 8, height = 4, units = "in", dpi = "retina")
