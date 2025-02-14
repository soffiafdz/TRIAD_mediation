#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
library(stringr)
library(lubridate)
library(glue)
library(ggplot2)

ReDoPlots   <- FALSE

### INPUT
fpaths      <- here("data", c("/derivatives/vols_hcvc_adni.csv",
                              "/rds/adni_dxs_imputed.rds",
                              "/rds/adni_age_calculated.rds"))

## HcVc VOLUMES
if (!file.exists(fpaths[1])) {
  sprintf("File: %s is required but could not be found.", fpath) |> stop()
}

hcvc.dt     <- fread(fpaths[1])

## Clinical DIAGNOSES
if (file.exists(fpaths[2])) {
  dx.dt     <- read_rds(fpaths[2])
} else {
  here("code/impute_dx.R") |> source()
}

## AGE
## Calculated age using the date of MRI sessions
if (file.exists(fpaths[3])) {
  age.dt    <- read_rds(fpaths[3])
} else {
  here("code/calculate_age_adni.R") |> source()
}

### CLEAN and MERGE
## Extract PTID and EXAMDATE from the filename of HcVc volumes
hcvc.dt[, PTID := str_extract(ID, "\\d{3}_S_\\d{4}")]
hcvc.dt[, EXAMDATE := ymd(str_extract(ID, "\\d{4}(-\\d{2}){2}|\\d{8}"))]
fnames.dt   <- hcvc.dt[, .(PTID, EXAMDATE, FNAME = ID)] |>
  setkey(PTID, EXAMDATE)
hcvc.dt[, ID := NULL]
setkey(hcvc.dt, PTID, EXAMDATE)

hcvc.dt     <- age.dt[, -"VISCODE"][hcvc.dt]
hcvc.dt     <- dx.dt[, -"VISCODE"][hcvc.dt]
rm(dx.dt, age.dt)

### IQR
## function
IQR.outliers <- function(x, side = "all") {
  if(any(is.na(x))) stop("Missing values found.")
  if(!is.numeric(x)) stop("Not numeric.")
  opt <- match.arg(side, choices = c("low", "high", "all"))
  Q1  <- quantile(x, 0.25)
  Q3  <- quantile(x, 0.75)
  IQR <- IQR(x)
  lb  <- (Q1 - (1.5 * IQR))
  hb  <- (Q3 + (1.5 * IQR))
  if (opt == "low") {
    return(x[x < lb])
  } else if (opt == "high") {
    return(x[x > hb])
  } else {
    return(c(x[x < lb], x[x > hb]))
  }
}

## Ignore full volumes; focus on L/R
hcvc.dt[, c("HC", "CSF") := NULL]

## Melt dt -> dt.l
hcvc.dt.l   <- melt(hcvc.dt, measure = patterns("^L|R"))
rois        <- c("LHC", "RHC", "LCSF", "RCSF")
# Voxel -> CC
hcvc.dt.l[, value := value / 1000]
for (roi in rois) {
  out_low   <- hcvc.dt.l[variable == roi, IQR.outliers(value, "low")]
  out_high  <- hcvc.dt.l[variable == roi, IQR.outliers(value, "high")]
  hcvc.dt.l[variable == roi & value %in% out_low, OUTLIER := glue("{roi}-")]
  hcvc.dt.l[variable == roi & value %in% out_high, OUTLIER := glue("{roi}+")]
  for (dx in c("CN", "MCI", "Dementia")) {
    out_low   <- hcvc.dt.l[variable == roi & DX == dx,
                           IQR.outliers(value, "low")]
    out_high  <- hcvc.dt.l[variable == roi & DX == dx,
                           IQR.outliers(value, "high")]
    hcvc.dt.l[variable == roi & DX == dx & value %in% out_low,
              OUTLIER_dx := glue("{roi}-")]
    hcvc.dt.l[variable == roi & DX == dx & value %in% out_high,
              OUTLIER_dx := glue("{roi}+")]
  }
}

## Bring back to wide
hcvc.dt     <- hcvc.dt.l |>
  dcast(... ~ variable, value.var = c("value", "OUTLIER", "OUTLIER_dx"))

## Concatenate OUTLIER columns
regex       <- "\\w*(\\+|\\-)"

# General outliers
cols        <- paste("OUTLIER", c(rois), sep = "_")
hcvc.dt[, OUTLIER_ := do.call(paste, .SD), .SDcols = cols]
# FIRST ROUND: look for NO outliers rows
hcvc.dt[, OUTLIER := str_extract(OUTLIER_, regex)]
# SECOND ROUND: for non-NAs rows, extract ALL values
hcvc.dt[!is.na(OUTLIER),
        OUTLIER := apply(str_extract_all(OUTLIER_, regex, simplify = TRUE),
                         1, paste, collapse = ";")] |> suppressWarnings()
# Remove ",,+"
hcvc.dt[, OUTLIER := str_remove(OUTLIER, ";{2,}|;$")]

# By Dx
cols        <- paste("OUTLIER_dx", c(rois), sep = "_")
hcvc.dt[, OUTLIER_ := do.call(paste, .SD), .SDcols = cols]
# FIRST ROUND: look for NO outliers rows
hcvc.dt[, OUTLIER_dx := str_extract(OUTLIER_, regex)]
# SECOND ROUND: for non-NAs rows, extract ALL values
hcvc.dt[!is.na(OUTLIER_dx),
        OUTLIER_dx := apply(str_extract_all(OUTLIER_, regex, simplify = TRUE),
                         1, paste, collapse = ";")] |> suppressWarnings()
# Remove ",,+"
hcvc.dt[, OUTLIER_dx := str_remove(OUTLIER_dx, "\\;{2,}|\\;$")]

## RENAME/REMOVE columns
# Rename volume columns
cols        <- paste("value", rois, sep = "_")
hcvc.dt |> setnames(cols, str_sub(cols, 7))
# Remove OUTLIER_* columns
cols <- grep("_(L|R)|_$", names(hcvc.dt), value = TRUE)
hcvc.dt[, (cols) := NULL]
rm(cols)


### OUTPUT
# CSV with filenames for further QC
setkey(hcvc.dt, PTID, EXAMDATE)
outcsv      <- here("data/derivatives/adni_hcvc_outliers.csv")
fnames.dt[hcvc.dt
          ][!is.na(OUTLIER) | !is.na(OUTLIER_dx),
          .(FNAME, OUTLIER, OUTLIER_dx)] |>
  fwrite(outcsv)
rm(outcsv)


# Detailed data.table
outrds      <- here("data/rds/adni_hcvc_outliers.rds")
write_rds(hcvc.dt, outrds)
rm(outrds)

### PLOTS
## Further CLEANUP
hcvc.dt.l[, `:=`(SIDE = factor(str_sub(variable, end = 1),
                               levels = c("L", "R"),
                               labels = c("Left", "Right")),
                 ROI  = factor(str_sub(variable, 2),
                               levels = c("HC", "CSF"),
                               labels = c("Hippocampus", "Ventricle")),
                 OUTLIER = str_extract(OUTLIER, "\\-|\\+"),
                 OUTLIER_dx = str_extract(OUTLIER_dx, "\\-|\\+"))]
setnames(hcvc.dt.l, "value", "CC")
hcvc.dt.l[, variable := NULL]

## Palette
cbPalette    <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## HC vs HVR: Boxplots
f_plot1  <- here("plots/boxplot_hc-exploration")
fp1_png  <- paste0(f_plot1, ".png")
#fp1_tiff <- paste0(f_plot1, ".tiff")
rm(f_plot1)
if(!file.exists(fp1_png)
   #|| !file.exists(fp1_tiff)
   || ReDoPlots) {
  hcvc.dt.l |>
  ggplot(aes(x = ROI, y = CC, colour = SIDE)) +
    theme_classic(base_size = 12) +
    theme(
       text = element_text(size = 12),
       axis.text.y = element_text(size = 10),
       axis.text.x = element_text(size = 10),
       legend.position = "bottom") +
    geom_boxplot(outlier.shape = 21,
                 outlier.colour = "grey",
                 outlier.alpha = .7) +
    scale_colour_manual(values = cbPalette[c(2:3)]) +
    labs(title = "Exploration of segmentations",
         y = "Volume (CC)", x = "Region", colour = "Side")

  if(!file.exists(fp1_png) || ReDoPlots){
    ggsave(fp1_png, width = 5, height = 5, units = "in", dpi = 600)
  }
  rm(fp1_png)

  #if(!file.exists(fp1_tiff) || ReDoPlots){
    #ggsave(fp1_tiff, width = 5, height = 5, units = "in",
           #device = "tiff", dpi = 600)
  #}
  #rm(fp1_tiff)
}

## HC vs HVR: Boxplots by Dx
f_plot2  <- here("plots/boxplot_hc-exploration-dx")
fp2_png  <- paste0(f_plot2, ".png")
#fp2_tiff <- paste0(f_plot2, ".tiff")
rm(f_plot2)
if(!file.exists(fp2_png)
   #|| !file.exists(fp2_tiff)
   || ReDoPlots) {
  hcvc.dt.l[!is.na(DX),
            .(CC, SIDE,
              ROI = factor(ROI, labels = c("Hipp", "Vent")),
              DX = str_remove(DX, "\\?"))] |>
  ggplot(aes(x = ROI, y = CC, colour = SIDE)) +
    theme_classic(base_size = 12) +
    theme(
       text = element_text(size = 12),
       axis.text.y = element_text(size = 10),
       axis.text.x = element_text(size = 10),
       legend.position = "bottom") +
    geom_boxplot(outlier.shape = 21,
                 outlier.colour = "grey",
                 outlier.alpha = .7) +
    facet_grid(cols = vars(DX)) +
    scale_colour_manual(values = cbPalette[c(2:3)]) +
    labs(title = "Exploration of segmentations by Diagnosis",
         y = "Volume (CC)", x = "Region", colour = "Side")

  if(!file.exists(fp2_png) || ReDoPlots){
    ggsave(fp2_png, width = 7, height = 5, units = "in", dpi = 600)
  }
  rm(fp2_png)

  #if(!file.exists(fp1_tiff) || ReDoPlots){
    #ggsave(fp1_tiff, width = 5, height = 5, units = "in",
           #device = "tiff", dpi = 600)
  #}
  #rm(fp2_tiff)
}

### HC vs HVR: Points
f_plot3  <- here("plots/pointplot_hc-exploration")
fp3_png  <- paste0(f_plot3, ".png")
#fp3_tiff <- paste0(f_plot3, ".tiff")
rm(f_plot3)
if(!file.exists(fp3_png)
   #|| !file.exists(fp3_tiff)
   || ReDoPlots) {
  dt <- copy(hcvc.dt.l)
  dt[, DX := str_remove(DX, "\\?")]
  dt[is.na(AGE), `:=`(AGE = 45, MISSING_age = TRUE)]
  lims.dt <- dt[, .(lim_low = quantile(CC, 0.25) - (1.5 * IQR(CC)),
                    lim_high = quantile(CC, 0.75) + (1.5 * IQR(CC))), ROI]
  lims_dx.dt <- dt[!is.na(DX),
                   .(lim_low = quantile(CC, 0.25) - (1.5 * IQR(CC)),
                     lim_high = quantile(CC, 0.75) + (1.5 * IQR(CC))),
                   .(ROI, DX)]
  dt |>
  ggplot(aes(x = AGE, y = CC, colour = DX)) +
    theme_classic(base_size = 12) +
    theme(
       text = element_text(size = 12),
       axis.text.y = element_text(size = 10),
       axis.text.x = element_text(size = 10),
       legend.position = "bottom") +
    geom_point(aes(shape = SIDE), size = 1, fill = "transparent", alpha = .7) +
    geom_hline(data = lims.dt, aes(yintercept = lim_low),
               colour = "grey", alpha = .7, linetype = "dashed") +
    geom_hline(data = lims.dt, aes(yintercept = lim_high),
               colour = "grey", alpha = .7, linetype = "dashed") +
    geom_hline(data = lims_dx.dt, aes(yintercept = lim_low, colour = DX),
               alpha = .4, linetype = "dashed") +
    geom_hline(data = lims_dx.dt, aes(yintercept = lim_high, colour = DX),
               alpha = .4, linetype = "dashed") +
    #geom_label_repel(dt[!is.na(SEX)][order(VAL), .SD[.N/2], .(DX, SEX, HC_msr)],
                     #aes(y = kappa,
                         #label = glue("{round(mean, 2)}\n({round(sd, 2)})")),
                     #size = 3, box.padding = 1.5, alpha = 0.75) +
    scale_fill_manual(values = c("white", "black"),
                      labels = c("Missing Age", "Calc. Age")) +
    scale_shape_manual(values = 21:23) +
    scale_colour_manual(values = cbPalette[c(2:4)]) +
    facet_wrap(facets = vars(ROI), nrow = 1, scales = "free") +
    labs(title = "Longitudinal exploration of segmentations for outliers",
         subtitle = "dashed lines represent +/- 1.5 * IQR; general & by Dx.",
         y = "Volume (CC)", x = "Age (years)",
         colour = "Dx", shape = "Side")

  if(!file.exists(fp3_png) || ReDoPlots){
    ggsave(fp3_png, width = 15, height = 15, units = "in", dpi = 200)
  }

  #if(!file.exists(fp1_tiff) || ReDoPlots){
    #ggsave(fp1_tiff, width = 5, height = 5, units = "in",
           #device = "tiff", dpi = 600)
  #}
}
