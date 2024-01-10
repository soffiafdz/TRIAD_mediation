#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lubridate)
library(readr)
library(gtsummary)

## Read/Parse CSV files
fpaths      <- here("data", c("demographics.csv",
                              "neuropsych_eval.csv",
                              "pet_biomarkers.csv"))

if (any(!file.exists(fpaths))) {
  here("code/parse_csv_data.R") |> source()
} else {
  demog.dt    <- fread(fpaths[1])
  neuropsy.dt <- fread(fpaths[2])
  pet.dt      <- fread(fpaths[3])
}

rm(fpaths)

## Read/source HVR
#vols.rds    <- here("data/rds", paste("hcv_hvr_adj-", c("all", "old"), ".rds"))
# Pick between controlling with old subset or old + youth
#vols.rds    <- here("data/rds/hcv_hvr_adj-all.rds") # Controlled by Old + Youth
vols.rds    <- here("data/rds/hcv_hvr_adj-old.rds") # Controlled by just old
if (file.exists(vols.rds)) {
  vols.dt <- read_rds(vols.rds)
} else {
  here("code/calc_hvr.R") |> source()
  #vols.dt <- vols_all.dt # Youth & Old
  vols.dt <- vols_old.dt # Just Old
}
rm(vols.rds)

## Assign VISIT labels to vols.dt
#VM00
vm00.dt     <- vols.dt[, .SD[which.min(SCANDATE)], PTID
                       ][, .(PTID, SCANDATE, VISIT = "VM00")]
#VM06
vm06.dt     <- vols.dt[!vm00.dt, on = .(PTID, SCANDATE)
                       ][PTID %in% demog.dt[VISIT == "VM06", PTID],
                       .SD[which.min(SCANDATE)], PTID
                       ][, .(PTID, SCANDATE, VISIT = "VM06")]
#VM12
vm12.dt     <- vols.dt[!vm00.dt, on = .(PTID, SCANDATE)
                       ][!vm06.dt, on = .(PTID, SCANDATE)
                       ][PTID %in% demog.dt[VISIT == "VM12", PTID],
                       .SD[which.min(SCANDATE)], PTID
                       ][, .(PTID, SCANDATE, VISIT = "VM12")]
#VM24
vm24.dt     <- vols.dt[!vm00.dt, on = .(PTID, SCANDATE)
                       ][!vm06.dt, on = .(PTID, SCANDATE)
                       ][!vm12.dt, on = .(PTID, SCANDATE)
                       ][PTID %in% demog.dt[VISIT == "VM24", PTID],
                       .SD[which.min(SCANDATE)], PTID
                       ][, .(PTID, SCANDATE, VISIT = "VM24")]
#VM36
vm36.dt     <- vols.dt[!vm00.dt, on = .(PTID, SCANDATE)
                       ][!vm06.dt, on = .(PTID, SCANDATE)
                       ][!vm12.dt, on = .(PTID, SCANDATE)
                       ][!vm24.dt, on = .(PTID, SCANDATE)
                       ][PTID %in% demog.dt[VISIT == "VM36", PTID],
                       .SD[which.min(SCANDATE)], PTID
                       ][, .(PTID, SCANDATE, VISIT = "VM36")]

visits.dt   <- rbindlist(list(vm00.dt, vm06.dt, vm12.dt, vm24.dt, vm36.dt))

vols.dt     <- visits.dt[vols.dt, on = .(PTID, SCANDATE)]
rm(visits.dt, vm00.dt, vm06.dt, vm12.dt, vm24.dt, vm36.dt)

# Missing subject is MRT62, VM24
vols.dt[is.na(VISIT), `:=`(PTID = "MRT62", VISIT = "VM24")]

# Imaging
imag.dt     <- pet.dt[vols.dt, on = .(PTID, VISIT)]

# COVARS
# Repeated MRT63 with NAs
covars.dt   <- neuropsy.dt[!(PTID == "MRT63" & is.na(EVALDATE))
                           ][imag.dt, on = .(PTID, VISIT)]

# All baseline data
triad.dt    <- demog.dt[covars.dt, on = .(PTID, VISIT)]
#rm(vols.dt, pet.dt, imag.dt, neuropsy.dt, covars.dt, demog.dt)

## Data cleaning
# Sex
triad.dt[, SEX := factor(SEX, labels = c("Female", "Male"))]

# Time differences
triad.dt[, AGE_scan := interval(ymd(DOB), ymd(SCANDATE)) / years(1)]
triad.dt[, EVAL_delay := ymd(EVALDATE) - ymd(SCANDATE)]
# MRT62 lacking SCANDATE use EVALDATE
triad.dt[is.na(AGE_scan),
         `:=`(AGE_scan = interval(ymd(DOB), ymd(EVALDATE)) / years(1),
              EVAL_delay = 0)]

# Clean DX
triad.dt[DX == "Unknown", DX := NA]
triad.dt[, DX_clean := DX]
triad.dt[DX %in% c("Atypical Dementia", "FTD", NA), DX_clean := "Other"]

# Remove useless columns
triad.dt[, c("DOB", "SCANDATE", "EVALDATE") := NULL]

# Export
write_rds(triad.dt, here("data/rds/triad.rds"))

# TODO: Decide if remove NAs
# Must have: Imaging
triad.dt    <- triad.dt[!is.na(AMYLOID) & !is.na(TAU_braak1) & !is.na(HVR_l)]

# Must have: Neuropsy
#triad.dt    <- triad.dt[!is.na(RAVLT_raw) & !is.na(MOCA_score)]

# No of evals
sessn       <- triad.dt[, .(SESS = as.character(.N)), PTID]
triad_bl.dt <- triad.dt[VISIT == "VM00"][sessn, on = "PTID"]
rm(sessn)

## Table1
fname       <- here("data/derivatives/table1.docx")
if (!file.exists(fname)) {
  triad_bl.dt[!is.na(APOE_n) & !is.na(MOCA_score) & DX %in% c("CN", "MCI"),
              .(DX_clean, SEX, AGE_scan, EDUC, APOE = factor(APOE_n),
                MOCA_score, SESS,
                #RAVLT_intro, RAVLT_raw, RAVLT_rep,

                AMYLOID, TAU_sum = (TAU_braak1 + TAU_braak2 + TAU_braak3 +
                                    TAU_braak4 + TAU_braak5 + TAU_braak6),
                #HCv_l, HCv_r, HVR_l, HVR_r)] |>
                HVR = (HVR_l + HVR_r) /2)] |>
    tbl_summary(by = DX_clean,
                label = list(SEX ~ "Sex",
                             AGE_scan ~ "Age (years)",
                             EDUC ~ "Education (years)",
                             APOE ~ "APOE4 alleles",
                             MOCA_score ~ "MoCA score",
                             #RAVLT_raw ~ "RAVLT (raw score)",
                             #RAVLT_intro ~ "RAVLT (intro score)",
                             #RAVLT_rep ~ "RAVLT (rep score)",
                             SESS ~ "Number of visits",
                             AMYLOID ~ "Amyloid (PET)",
                             TAU_sum ~ "Tau (PET)"),
                             #HCv_l ~ "HC vol (left)",
                             #HCv_r ~ "HC vol (right)",
                             #HVR_l ~ "HVR (left)",
                             #HVR_r ~ "HVR (right)"),
                statistic = all_continuous() ~ "{mean} ({sd})",
                #missing_text = "Missing") |>
                missing = "no") |>
    modify_header(label ~ "**Variable**") |>
    #modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Clinical Label**") |>
    #add_n() |>
    add_p() |> as_flex_table() |>
    flextable::save_as_docx(path = fname)
}
