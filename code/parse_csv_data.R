#!/usr/bin/env Rscript

library(here)
library(data.table)

## Read files
# File names
fpaths <- c(
  "MCSA_demo_scan_transfer_ANTS_processing_new_20220206_anonym",
  "from_cecile_data-2023-02-10T16_35_19.435Z",
  "CEREBRA_volumetric_20231208"
) |>
  sprintf(fmt = "data/data_2023/%s.csv") |>
  here()

# Parsing
data.lst <- list(
  PET = fread(fpaths[1]),
  NEURO = fread(fpaths[2]),
  CEREBRA = fread(fpaths[3])
)
rm(fpaths)

## DTs
# Demographics
demog_cols  <- c(
  "FID", "visit", "DX_cat", "dob", "sex", "edu", "apoe_add", "WMH_wm", "MMSE"
)

demog.dt    <- data.lst[["PET"]][, ..demog_cols]

setnames(
  demog.dt,
  demog_cols,
  c("PTID", "VISIT", "DX", "DOB", "SEX", "EDUC", "APOE_n", "WMH", "MMSE")
)

fwrite(demog.dt, here("data/data_2023/demographics.csv"))

# PET biomarkers (Global/Braak Stage)
pet_cols    <- c(
  "FID", "visit", "NeoctxAZD_SUVR", # Amyloid: AZD4694
  # Tau BraakStages: MK-6240 ligand
  paste0(paste0("Braak", c(1:6, "Stage")), "_masked")
)

pet.dt      <- data.lst[["PET"]][, ..pet_cols]

setnames(
  pet.dt,
  pet_cols,
  c("PTID", "VISIT", "AMYLOID", paste0("TAU_braak", 1:6), "TAU_braak_stage")
)

fwrite(pet.dt, here("data/data_2023/pet_biomarkers.csv"))

# Amyloid & Tau by ROI (CEREBRA)
roi_cols    <- c(
  "FID", "visit", "label_id", "label", "side", "vol", "mk_pet", "nav_pet"
)

cerebra.dt  <- data.lst[["CEREBRA"]][, ..roi_cols]

setnames(
  cerebra.dt,
  roi_cols,
  c("PTID", "VISIT", "LABEL_id", "LABEL_name", "SIDE", "VOL", "TAU", "AMYLOID")
)

## Normalize SUVR values by Avg L/R cerebellar gray matter
cerebra.dt  <- cerebra.dt[
 "Cerebellum_Gray_Matter",
  on = "LABEL_name",
  .(
    TAU_cgm = mean(TAU),
    AMY_cgm = mean(AMYLOID)
  ),
  .(PTID, VISIT)
][
  cerebra.dt,
  on = .(PTID, VISIT),
  .(
    PTID, VISIT, LABEL_id, LABEL_name, SIDE, VOL,
    TAU_norm = TAU / TAU_cgm,
    AMYLOID_norm = AMYLOID / AMY_cgm
  )
]

fwrite(cerebra.dt, here("data/data_2023/pet_biomarkers_cerebra.csv"))

# Neuropsych
cog_cols    <- c(
  "FID",
  "FVIS",
  "MOCA_score",
  paste0(
    "Neuropsych_RAVLT_",
    c(
      "Date_taken",
      paste0("trial_B1_", c("intrusion", "raw", "repetition"), "_score")
    )
  )
)

neuropsy.dt <- data.lst[["NEURO"]][, ..cog_cols]

setnames(
  neuropsy.dt,
  cog_cols,
  c(
    "PTID",
    "VISIT",
    "MOCA_score",
    "EVALDATE",
    paste0("RAVLT_", c("intro", "raw", "rep"))
  )
)

fwrite(neuropsy.dt, here("data/data_2023/neuropsych_eval.csv"))

rm(demog_cols, pet_cols, roi_cols, cog_cols)
