#!/usr/bin/env Rscript

library(here)
library(data.table)

## Read files
# File names
fname1  <- here("data/MCSA_demo_scan_transfer_ANTS_processing_new_20220206_anonym.csv")
fname2  <- here("data/from_cecile_data-2023-02-10T16_35_19.435Z.csv")
fname3  <- here("data/CEREBRA_volumetric_20231208.csv")

# Parsing
dt1     <- fread(fname1)
dt2     <- fread(fname2)
dt3     <- fread(fname3)
rm(fname1, fname2, fname3)

## DTs
# Demographics
demog_cols  <- c("FID", "visit", "DX_cat",
                 "dob", "sex", "edu", "apoe_add", "WMH_wm", "MMSE")
demog.dt    <- dt1[, ..demog_cols]
setnames(demog.dt, demog_cols,
         c("PTID", "VISIT", "DX", "DOB", "SEX", "EDUC", "APOE_n", "WMH",
           "MMSE"))
fwrite(demog.dt, here("data/demographics.csv"))
rm(demog_cols)

# PET biomarkers (Global/Braak Stage)
pet_cols    <- c("FID", "visit", "NeoctxAZD_SUVR", # Amyloid: AZD4694
                 # Tau BraakStages: MK-6240 ligand
                 paste0(paste0("Braak", c(1:6, "Stage")), "_masked"))
pet.dt      <- dt1[, ..pet_cols]
setnames(pet.dt, pet_cols, c("PTID", "VISIT", "AMYLOID",
                             paste0("TAU_braak", 1:6), "TAU_braak_stage"))
fwrite(pet.dt, here("data/pet_biomarkers.csv"))
rm(pet_cols, dt1)

# Amyloid & Tau by ROI (CEREBRA)
roi_cols    <- c("FID", "visit", "label_id", "label", "side",
                 "vol", "mk_pet", "nav_pet")
cerebra.dt  <- dt3[, ..roi_cols]
setnames(cerebra.dt, roi_cols, c("PTID", "VISIT", "LABEL_id", "LABEL_name",
                                 "SIDE", "VOL", "TAU", "AMYLOID"))
fwrite(cerebra.dt, here("data/pet_biomarkers_cerebra.csv"))
rm(roi_cols, dt3)

# Neuropsych
cog_cols    <- c("FID", "FVIS", "MOCA_score",
                 paste0("Neuropsych_RAVLT_", c("Date_taken",
                                               paste0("trial_B1_",
                                                      c("intrusion", "raw",
                                                        "repetition"),
                                                      "_score"))))
neuropsy.dt <- dt2[, ..cog_cols]
setnames(neuropsy.dt, cog_cols, c("PTID", "VISIT", "MOCA_score", "EVALDATE",
                                  paste0("RAVLT_", c("intro", "raw", "rep"))))
fwrite(neuropsy.dt, here("data/neuropsych_eval.csv"))
rm(cog_cols, dt2)
