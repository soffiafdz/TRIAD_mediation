#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
library(lme4)
library(ggplot2)
library(patchwork)
library(ADNIMERGE)

ReDoPlots   <- TRUE
USE_IMPUTED <- TRUE

### INPUT
fpaths      <- here("data/rds", c("adni_hc-hvr.rds",
                                  "adni_dxs_imputed.rds",
                                  "adni_age_calculated.rds",
                                  "ucb_pet-amy.rds",
                                  "ucb_pet-tau.rds",
                                  "../adni_wmh-vol.csv",
                                  "adni_cog-latent.rds",
                                  "adni_biomarkers-latent.rds"))

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

## PET
if (all(file.exists(fpaths[4:5]))) {
  ucb_amy.dt <- readRDS(fpaths[4])
  ucb_tau.dt <- readRDS(fpaths[5])
} else {
  here("code/parse_pet_adni.R") |> source()
}

## WMH
if (!file.exists(fpaths[6])) {
  sprintf("File: %s is required but could not be found.", fpaths[6]) |> stop()
}
wmh.dt      <- fread(fpaths[6])

## Cognition
# Latent variable obtained from ADASQ4, MMSE, CDRSB and validated with CFA
if (file.exists(fpaths[7])) {
  cog.dt    <- read_rds(fpaths[7])
} else {
  here("code/cfa-cog_adni.R") |> source()
}

## Pathology
# Latent variables obtained from UCBerkeley Amy, Tau PET and WMH vols
if (file.exists(fpaths[8])) {
  path.dt   <- read_rds(fpaths[8])
} else {
  here("code/cfa-biomarkers_adni.R") |> source()
}
rm(fpaths)


### Data PROCESSING
## Use imputed Dx?
Dxs         <- c("CN", "Dementia")
if (USE_IMPUTED) dx.dt[, DX := stringr::str_remove(DX, "\\?")]
sel_dx.dt <- dx.dt[, DX, .(PTID, EXAMDATE)]

## PTGENDER & APOE4
data(adnimerge)
setDT(adnimerge)
tiv.dt      <- adnimerge[!duplicated(PTID), .(SEX = PTGENDER, APOE4), PTID]

## HC-HVR
# MERGE and average by side
DT          <- tiv.dt[age.dt, on = "PTID"
                      ][sel_dx.dt, on = .(PTID, EXAMDATE)
                      ][hc_hvr.dt, on = .(PTID, EXAMDATE)
                      ][, .(HCv = mean(HCvol_adj), HVR = mean(HVR)),
                      .(PTID, EXAMDATE, DX, AGE, SEX, APOE4)
                      ][wmh.dt, on = .(PTID, EXAMDATE)
                      ][ucb_amy.dt[, c(1:2,4,7)],
                      on = .(PTID, EXAMDATE = DATE_MRI)
                      ][ucb_tau.dt[, c(1:2,4,7)],
                      on = .(PTID, EXAMDATE = DATE_MRI)
                      ][path.dt, on = .(PTID, EXAMDATE)] |>
na.omit() |>
setkey(PTID)
rm(tiv.dt, hc_hvr.dt, wmh.dt, ucb_amy.dt, ucb_tau.dt, path.dt)

## Data CLEANING
# Do this functionally because the sample will be changing depending on the
# Covariates of interest
DTclean     <- function(DT, time = TRUE, scalevars = NULL, ordervars = NULL,
                        scale_baseline = TRUE) {
  if (time) {
    DTw <- copy(DT)
    DTw[, c("AGE.bl", "TIME") := NULL] |> suppressWarnings()
    DTw[order(AGE), ID := rowid(PTID)]
    setcolorder(DTw, "ID", after = "PTID")
    DTw <- DTw[DTw[.(1), on = "ID", .(AGE.bl = AGE), "PTID"]]
    DTw[, TIME := AGE - AGE.bl]
    setcolorder(DTw, c("AGE.bl", "TIME"), before = "AGE")
  }
  if (!is.null(ordervars)) {
    DTw[, (ordervars) := lapply(.SD, as.ordered), .SDcols = ordervars]
  }
  if (!is.null(scalevars)) {
    if (is.integer(scalevars)) scalevars <- names(DTw)[scalevars]
    if (scale_baseline) {
      DTp <- rbind(DTw[.(1), on = "ID", lapply(.SD, mean), .SDcols = scalevars],
                   DTw[.(1), on = "ID", lapply(.SD, sd), .SDcols = scalevars])
      DTw[, (paste0(scalevars, "_scl")) :=
          lapply(scalevars, \(x) {(get(x) - DTp[[x]][1]) / DTp[[x]][2]})]
    } else {
      DTw[, (paste0(scalevars, "_scl")) := lapply(.SD, scale),
          .SDcols = scalevars]
    }
  }
  return(DTw)
}

DTunscale <- function(DT, scalevars, origs,
                      scale_baseline = TRUE, replace = TRUE, newnames = NULL) {
  DTw <- copy(DT)
  if (scale_baseline) {
    DTp <- rbind(DTw[.(1), on = "ID", lapply(.SD, mean), .SDcols = origs],
                 DTw[.(1), on = "ID", lapply(.SD, sd), .SDcols = origs])
  } else {
    DTp <- rbind(DTw[, lapply(.SD, mean), .SDcols = origs],
                 DTw[, lapply(.SD, sd), .SDcols = origs])
  }
  setnames(DTp, scalevars)
  if (replace) {
    DTw[, (scalevars) :=
        lapply(scalevars, \(x) {(get(x) * DTp[[x]][2]) + DTp[[x]][1]})]
  } else {
    newnames <- if (!is.null(newnames)) {
      newnames
    } else if (any(stringr::str_detect(scalevars, "scl"))) {
      sub("scl", "uscl", scalevars)
    } else {
      paste0(scalevars, "_uscl")
    }
    DTw[, (newnames) :=
        lapply(scalevars, \(x) {(get(x) * DTp[[x]][2]) + DTp[[x]][1]})]
  }
}

setnames(DT, c("CENTILOIDS", "META_TEMPORAL_SUVR"), c("AMY", "TAU"))
#setcolorder(DT, c(1:8, 10, 12:15))
#DT          <- DTclean(DT, scalevars = colnames(DT)[6:13], ordervars = "APOE4")

### PLOTS
### lme4
#f <- paste("%s ~ TIME * TAU_scl + AMY_scl",
           #"+ SEX + AGE.bl + APOE4 + I(TIME^2) + (1|PTID)") |>
#sprintf(c("HCv_scl", "HVR_scl")) |>
#lapply(as.formula)
#hcv.mod <- lmer(f[[1]], DT, REML = F,
                #control = lmerControl(optimizer = "bobyqa"))
#hvr.mod <- lmer(f[[2]], DT, REML = F,
                #control = lmerControl(optimizer = "bobyqa"))

#DT[, HCvp := predict(hcv.mod)]
#DT[, HVRp := predict(hvr.mod)]
#DT      <- DTunscale(DT, scalevars = c("HCvp", "HVRp"), orig = c("HCv", "HVR"))
##DT      <- melt(DT, measure = patterns(OBS = "scl$", PRED = "p$"),
                ##variable.name = "HC")
##DT[, HC := factor(HC, labels = c("Hipp. Vol.", "HVR"))]
##DT[, DX := factor(DX, levels = c("CN", "MCI", "Dementia"),
                  ##labels = c("CN", "MCI", "AD"))]

#### Cognition
#DT_c        <- age.dt[cog.dt, on = .(PTID, VISCODE), .(PTID, AGE, COG_latent)
                      #][DT, on = .(PTID, AGE)] |>
#na.omit() |>
#setkey(PTID)

#DTc <- copy(DT_c)
##setcolorder(DT_c, c(1:2,4:6))
#DT_c[, COG_orig := COG_latent]
#DT_c[, COG_latent := COG_latent + abs(min(COG_latent))]
##DT_c        <- DTclean(DT_c, ordervars = "APOE4", scalevars = names(DT_c)[6:14])

### MODELS
##f <- paste("COG_latent ~ %s * TIME * APOE4 + I(TIME^2)",
           ##"* APOE4 + AGE.bl + SEX + (TIME|PTID)") |>
##sprintf(c("HCv_scl", "HVR_scl")) |>
##lapply(as.formula)
##hcv_c.mod <- lmer(f[[1]], DT_c, REML = F,
                  ##control = lmerControl(optimizer = "bobyqa"))

##hvr_c.mod <- lmer(f[[2]], DT_c, REML = F,
                  ##control = lmerControl(optimizer = "bobyqa"))

##DT_c[, HCv_cog := predict(hcv_c.mod)]
##DT_c[, HVR_cog := predict(hvr_c.mod)]

##DT_c    <- melt(DT_c, measure = patterns(SCALED = "scl$", PRED = "cog$"),
                ##variable.name = "HCvar")

##DT_c[, HCvar := factor(HCvar, labels = c("Volume", "HVR"))]
##DT_c <- DT[, DX, .(PTID, AGE)][DT_c, on = .(PTID, AGE)] |> unique()


#### Palette
##cbPalette    <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                 ##"#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### HC vs HVR
##f_plot1  <- here("plots/adni_trajectory_hcv-hvr_mlm")
##fp1_png  <- paste0(f_plot1, ".png")
###fp1_tiff <- paste0(f_plot1, ".tiff")
##rm(f_plot1)
##if(!file.exists(fp1_png)
   ###|| !file.exists(fp3_tiff)
   ##|| ReDoPlots) {
  ##ggplot(DT, aes(x = AGE, y = OBS, colour = DX, fill = DX)) +
    ##theme_classic(base_size = 12) +
    ##theme(
       ##text = element_text(size = 12),
       ##axis.text.y = element_text(size = 10),
       ##axis.text.x = element_text(size = 10),
       ##legend.position = "bottom") +
    ###geom_point(data = DT[, .SD[which.max(ID)], PTID][ID == 1],
               ###shape = 21, size = .3, fill = "transparent", alpha = .5) +
    ##geom_line(aes(group = PTID), size = .3, linewidth = .3, alpha = .3) +
    ##geom_smooth(aes(y = PRED), size = .75, linewidth = .5, alpha = .2) +
    ##scale_colour_manual(values = cbPalette[c(2,8,6)]) +
    ##scale_fill_manual(values = cbPalette[c(2,8,6)]) +
    ##facet_wrap(facets = vars(HC), nrow = 1) +
    ##labs(title = "Neurodegeneration trajectories:",
         ##subtitle = "Mixed-effect models: Hippocampal volume & HVR.",
         ##y = "SD from CN mean", x = "Age (years)", colour = "Dx", fill = "Dx")

  ##if(!file.exists(fp1_png) || ReDoPlots){
    ##ggsave(fp1_png, width = 6, height = 6, units = "in", dpi = 600)
  ##}

  ###if(!file.exists(fp1_tiff) || ReDoPlots){
    ###ggsave(fp1_tiff, width = 5, height = 5, units = "in",
           ###device = "tiff", dpi = 600)
  ###}
##}

##f_plot2  <- here("plots/adni_trajectory_cog-lat")
##fp2_png  <- paste0(f_plot2, ".png")
###fp2_tiff <- paste0(f_plot2, ".tiff")
##rm(f_plot2)
##if (T) {
###if(!file.exists(fp2_png)
   ####|| !file.exists(fp3_tiff)
   ###|| ReDoPlots) {
  ##p1 <- ggplot(DT_c, aes(x = AGE, y = COG_latent, colour = DX)) +
    ##theme_classic(base_size = 12) +
    ##theme(
       ##text = element_text(size = 12),
       ##axis.text.y = element_text(size = 10),
       ##axis.text.x = element_text(size = 10),
       ##legend.position = "none") +
    ##geom_line(aes(group = PTID), size = .3, linewidth = .3, alpha = .2) +
    ##geom_hline(yintercept = 0, lty = "dashed", colour = cbPalette[1]) +
    ##scale_colour_manual(values = cbPalette[c(2,8,6)]) +
    ##labs(subtitle = "Individual cognitive decline trajectories",
         ##y = "Decline (observed)", x = "Age (years)")

  ##p2 <- DT_c[AGE > 57 & AGE < 93]["HVR", on = "HCvar"] |>
  ##ggplot(aes(x = AGE, y = COG_latent)) +
    ##theme_classic(base_size = 12) +
    ##theme(
       ##text = element_text(size = 12),
       ##axis.text.y = element_text(size = 10),
       ##axis.text.x = element_text(size = 10),
       ##legend.position = "right") +
    ###geom_smooth(colour = cbPalette[1], size = .5, linewidth = .5, se = T, method = "lm") +
    ##scale_fill_manual(values = cbPalette[c(2,8,6)]) +
    ##scale_colour_manual(values = cbPalette[c(2,8,6)]) +
    ##geom_hline(yintercept = 0, lty = "dashed", colour = cbPalette[1]) +
    ##geom_smooth(aes(colour = DX, fill = DX), size = .3,
                ##linewidth = .5, se = T, alpha = .3) +
    ##geom_smooth(colour = cbPalette[1], fill = cbPalette[1], size = .5,
                ##linewidth = .7, se = T, alpha = .3) +
    ##labs(subtitle = "Group cognitive decline trajectories",
         ##y = "Decline (predicted)", x = "Age (years)",
         ##colour = "Dx", fill = "Dx")

  ##p3 <- DT_c |>
    ##ggplot(aes(x = COG_latent, y = SCALED, colour = HCvar, fill = HCvar)) +
    ##theme_classic(base_size = 12) +
    ##theme(
       ##text = element_text(size = 12),
       ##axis.text.y = element_text(size = 10),
       ##axis.text.x = element_text(size = 10),
       ##legend.position = "right") +
    ##geom_point(size = .5, shape = 21, fill = "transparent", alpha = .15) +
    ##geom_smooth(size = .3, linewidth = .3, alpha = .07) +
    ##scale_shape_manual(values = 21:23) +
    ##scale_colour_manual(values = cbPalette[c(4,7)]) +
    ##scale_fill_manual(values = cbPalette[c(4,7)]) +
    ##labs(subtitle = "Cognitive Decline & Hippocampus",
         ##y = "HC (scaled)", x = "Cognitive decline",
         ##colour = "HC measure", fill = "HC measure")


  ###DT_c[, RESID := COG_latent - PRED]
  ###p3 <- ggplot(DT_c, aes(x = AGE, y = RESID, colour = HCvar)) +
    ###theme_classic(base_size = 12) +
    ###theme(
       ###text = element_text(size = 12),
       ###axis.text.y = element_text(size = 10),
       ###axis.text.x = element_text(size = 10),
       ###legend.position = "bottom") +
    ###geom_smooth(se = F, linewidth = .05, alpha = .5, method = "lm") +
    ####geom_smooth(se = F, linewidth = .05, alpha = .5) +
    ###scale_colour_manual(values = cbPalette[c(4,7)]) +
    ###labs(subtitle = "Mixed-effect model residuals",
         ###y = "Residuals", x = "Age (years)",
         ###colour = "Regression")

##p <- p1 / (p2 + p3) +
  ##plot_annotation(title = "Cognitive Decline Trajectories",
                  ##subtitle = paste("Mixed-effect regressions:",
                                   ##"Hippocampal volume & HVR;",
                                   ##"random intercepts and slopes"),
                  ##caption = paste("Cognitive decline was measured",
                                  ##"as a latent variable obtained from:",
                                  ##"ADASQ4, MMSE & CDRSB")) +
  ##plot_layout(guides = "collect")

  ##if(!file.exists(fp2_png) || ReDoPlots){
    ##ggsave(fp2_png, plot = p, width = 12, height = 8, units = "in", dpi = 600)
  ##}

  ###if(!file.exists(fp2_tiff) || ReDoPlots){
    ###ggsave(fp2_tiff, width = 5, height = 5, units = "in",
           ###device = "tiff", dpi = 600)
  ###}
##}
