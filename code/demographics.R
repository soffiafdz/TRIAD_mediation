#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lubridate)
library(readr)
library(gtsummary)

## Recreate tables
redo_tables <- TRUE

## Read/Parse CSV files
fpaths      <- here("data/rds", c("covars.rds", "raket_eds.rds",
                                  "wmh_vols.rds",
                                  "pet_cerebra.rds"))

if (!file.exists(fpaths[4]))   here("code/parse_pet.R") |> source()
if (any(!file.exists(fpaths)))  here("code/parse_csv_data.R") |> source()

covars.dt   <- read_rds(fpaths[1]) |> setkey(PTID, VISIT)
raket.dt    <- read_rds(fpaths[2]) |> setkey(PTID, VISIT)
wmh.dt      <- read_rds(fpaths[3]) |> setkey(PTID, VISIT)
pet.dt      <- read_rds(fpaths[4]) |> setkey(PTID, VISIT)

rm(fpaths)

## Read/source HVR
# Pick between controlling with old subset or old + youth
#vols.rds    <- here("data/rds/hcv_hvr_adj-young-old.rds") # Controlled by Old + Youth
vols.rds    <- here("data/rds/hcv_hvr_adj-old.rds") # Controlled by just old

if (file.exists(vols.rds)) {
  vols.dt <- read_rds(vols.rds)
} else {
  here("code/calc_hvr.R") |> source()
  if (grepl("young", vols.rds)) {
    vols.dt <- vols_all.dt # Youth & Old
  } else {
    vols.dt <- vols_old.dt # Just Old
  }
}
rm(vols.rds)

# Filter missing data
covars_subs.dt<- unique(covars.dt[, .(PTID, VISIT)])
pet_subs.dt   <- unique(pet.dt[, .(PTID, VISIT)])
vols_subs.dt  <- unique(vols.dt[, .(PTID, VISIT)])
wmh_subs.dt   <- unique(wmh.dt[, .(PTID, VISIT)])

raket_subs.dt <- unique(raket.dt[AB_bool == TRUE, .(PTID, VISIT)])

mri_subs.dt   <- wmh_subs.dt[vols_subs.dt, nomatch = 0]
imag_subs.dt  <- wmh_subs.dt[pet_subs.dt, nomatch = 0]

all_subs.dt   <- imag_subs.dt[raket_subs.dt, nomatch = 0]
write_rds(all_subs.dt, here("data/rds/incl_subs.rds"))

# Clean APOE4 data
covars.dt[APOE == "", APOE := NA]
covars.dt[APOE %in% c("E1/E2", "E3/E2"), APOE := "E2/E3"]
covars.dt[APOE == "E4/E2", APOE := "E2/E4"]
covars.dt[APOE == "E3//E3", APOE := "E3/E3"]
covars.dt[APOE == "E4/E3", APOE := "E3/E4"]

## Data merging
DT <-
  covars.dt[, .(PTID, VISIT, AGE, SEX, EDUC, APOE,
                CDR = as.numeric(CDR), MMSE = as.numeric(MMSE))
            ][raket.dt[AB_bool == TRUE, .(PTID, VISIT, RAKET_group)]
            ][pet.dt[, .(AMYLOID    = weighted.mean(SUVR_nav, VOLUME_nav),
                         TAU        = weighted.mean(SUVR_mk, VOLUME_mk)),
                     .(PTID, VISIT)]
            ][vols.dt[, .(PTID, VISIT,
                          HCv_mean = (HCv_l + HCv_r) / 2,
                          HVR_mean = (HVR_l + HVR_r) / 2)]
            ][wmh.dt
            ][all_subs.dt]


## Table1
fname       <- here("data/derivatives/table1_raket.docx")
if (!file.exists(fname) | redo_tables) {
  col_order <- c("RAKET_group", "AGE", "SEX", "EDUC", "APOE", "SESSN")

  DT[, .SD[which.min(AGE)], PTID
     ][DT[, .(SESSN = as.character(.N)), PTID]
     ][, ..col_order] |>
    tbl_summary(by = RAKET_group,
                label = list(SESSN ~ "Number of visits",
                             AGE ~ "Age (years)",
                             SEX ~ "Sex",
                             EDUC ~ "Education (years)",
                             APOE ~ "APOE4 alleles"),
                statistic = all_continuous() ~ "{mean} ({sd})",
                missing_text = "Missing") |>
                #missing = "no") |>
    modify_header(label ~ "**Variable**") |>
    modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Estimated disease stage**") |>
    add_p() |> as_flex_table() |>
    flextable::save_as_docx(path = fname)
}

fname       <- here("data/derivatives/table2_raket.docx")
if (!file.exists(fname) | redo_tables) {
  col_order <- c("RAKET_group", "MMSE", "WMH",
                 "AMYLOID", "TAU",
                 "HCv_mean", "HVR_mean")

  DT[, ..col_order] |>
    tbl_summary(by = RAKET_group,
                label = list(AMYLOID ~ "Amyloid (PET)",
                             TAU ~ "Tau (PET)",
                             HCv_mean ~ "HC vol",
                             HVR_mean ~ "HVR"
                             ),
                statistic = all_continuous() ~ "{mean} ({sd})",
                missing_text = "Missing") |>
                #missing = "no") |>
    modify_header(label ~ "**Variable**") |>
    modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Estimated disease stage**") |>
    add_p() |> as_flex_table() |>
    flextable::save_as_docx(path = fname)
}
