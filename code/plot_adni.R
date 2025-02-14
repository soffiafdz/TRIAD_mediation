#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
#library(stringr)
#library(lubridate)
#library(glue)
library(lme4)
library(ggplot2)
library(ADNIMERGE)

ReDoPlots   <- FALSE
USE_IMPUTED <- TRUE

### INPUT
fpaths      <- here("data/rds", c("adni_hc-hvr.rds",
                                  "adni_dxs_imputed.rds",
                                  "adni_age_calculated.rds"))

## HCvol & HVR
if (file.exists(fpaths[1])) {
  hc_hvr.dt <- read_rds(fpaths[1])
} else {
  here("code/calc_hvr_adni.R") |> source()
}

## Imputed Dxs from ADNIMERGE
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

rm(fpaths)


### Data PROCESSING
## Only CN or AD
## Use imputed Dx?
Dxs         <- c("CN", "Dementia")
if (USE_IMPUTED) dx.dt[, DX := stringr::str_remove(DX, "\\?")]
sel_dx.dt   <- dx.dt[DX %in% Dxs, DX, .(PTID, EXAMDATE)]
sel_dx.dt[, DX := factor(DX, labels = c("Cognitively healthy",
                                        "Dementia patients"))]
rm(USE_IMPUTED, dx.dt)

## PTGENDER
data(adnimerge)
setDT(adnimerge)
sex.dt      <- adnimerge[!duplicated(PTID), .(SEX = PTGENDER), PTID]

## MERGE and average by side
DT          <- sex.dt[age.dt, on = "PTID"
                      ][sel_dx.dt, on = .(PTID, EXAMDATE)
                      ][hc_hvr.dt, on = .(PTID, EXAMDATE)
                      ][!is.na(AGE)
                      ][!is.na(DX),
                      .(HCv = mean(HCvol_adj), HVR = mean(HVR)),
                      .(PTID, DX, AGE, SEX)]


### PLOTS
## lme4
mlm1 <- lmer(HCv ~ AGE * SEX * DX + (AGE|PTID), DT)
mlm2 <- lmer(HVR ~ AGE * SEX * DX + (AGE|PTID), DT)

## Palette
cbPalette    <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## HC vs HVR
f_plot1  <- here("plots/adni_trajectory_hcv_lme_rs")
fp1_png  <- paste0(f_plot1, ".png")
#fp1_tiff <- paste0(f_plot1, ".tiff")
rm(f_plot1)
if(!file.exists(fp1_png)
   #|| !file.exists(fp3_tiff)
   || ReDoPlots) {
  ggplot(DT, aes(x = AGE, y = HCv, colour = SEX)) +
    theme_classic(base_size = 12) +
    theme(
       text = element_text(size = 12),
       axis.text.y = element_text(size = 10),
       axis.text.x = element_text(size = 10),
       legend.position = "bottom") +
    geom_point(size = .5, fill = "transparent", alpha = .9) +
    geom_line(aes(y = predict(mlm1, re.form = NULL))) +
    #geom_line(aes(group = PTID), linewidth = .3, alpha = 5) +
    #geom_smooth(aes(fill = SEX), method = "lm", alpha = .3, linewidth = .3) +
    scale_colour_manual(values = cbPalette[c(8,6)]) +
    scale_fill_manual(values = cbPalette[c(8,6)]) +
    facet_wrap(facets = vars(DX), nrow = 1) +
    labs(title = "Neurodegeneration trajectories:",
         subtitle = "ADNI Cohort — Hippocampal volume.",
         y = "Volume (CC)", x = "Age (years)", colour = "Sex", fill = "Sex")

  if(!file.exists(fp1_png) || ReDoPlots){
    ggsave(fp1_png, width = 6, height = 6, units = "in", dpi = 600)
  }

  #if(!file.exists(fp1_tiff) || ReDoPlots){
    #ggsave(fp1_tiff, width = 5, height = 5, units = "in",
           #device = "tiff", dpi = 600)
  #}
}

f_plot2  <- here("plots/adni_trajectory_hvr_lme_rs")
fp2_png  <- paste0(f_plot2, ".png")
#fp2_tiff <- paste0(f_plot2, ".tiff")
rm(f_plot2)
if(!file.exists(fp2_png)
   #|| !file.exists(fp3_tiff)
   || ReDoPlots) {
  ggplot(DT, aes(x = AGE, y = HVR, colour = SEX)) +
    theme_classic(base_size = 12) +
    theme(
       text = element_text(size = 12),
       axis.text.y = element_text(size = 10),
       axis.text.x = element_text(size = 10),
       legend.position = "bottom") +
    geom_point(size = .5, fill = "transparent", alpha = .9) +
    geom_line(aes(y = predict(mlm2, re.form = NULL))) +
    #geom_line(aes(group = PTID), linewidth = .3, alpha = 5) +
    #geom_smooth(aes(fill = SEX), method = "lm", alpha = .3, linewidth = .3) +
    scale_colour_manual(values = cbPalette[c(8,6)]) +
    scale_fill_manual(values = cbPalette[c(8,6)]) +
    facet_wrap(facets = vars(DX), nrow = 1) +
    labs(title = "Neurodegeneration trajectories:",
         subtitle = "ADNI Cohort — Hippocampal-to-Ventricle Ratio (HVR).",
         y = "HVR (0-1)", x = "Age (years)", colour = "Sex", fill = "Sex")

  if(!file.exists(fp2_png) || ReDoPlots){
    ggsave(fp2_png, width = 6, height = 6, units = "in", dpi = 600)
  }

  #if(!file.exists(fp2_tiff) || ReDoPlots){
    #ggsave(fp2_tiff, width = 5, height = 5, units = "in",
           #device = "tiff", dpi = 600)
  #}
}
