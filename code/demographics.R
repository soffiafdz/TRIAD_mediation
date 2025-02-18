#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lubridate)
library(readr)
library(gtsummary)

## Recreate tables
REDOTABLES <- TRUE

## Read/Parse CSV files
fpaths      <- c(
  "demographics",
  "neuropsych_eval",
  "pet_biomarkers",
  "pet_biomarkers_cerebra"
) |>
  sprintf(fmt = "data/data_2023/%s.csv") |>
  here()

if (any(!file.exists(fpaths))) {
  here("code/parse_csv_data.R") |> source()
} else {
  demog.dt    <- fread(fpaths[1])
  neuropsy.dt <- fread(fpaths[2])
  pet.dt      <- fread(fpaths[3])
  cerebra.dt  <- fread(fpaths[4])
}

rm(fpaths)

## Read/source HVR
#vols.rds    <- here("data/rds", paste("hcv_hvr_adj-", c("all", "old"), ".rds"))
# Pick between controlling with old subset or old + youth
#vols.rds    <- here("data/rds/hcv_hvr_adj-all.rds") # Controlled by Old + Youth
vols.rds    <- here("data/rds/hcv_hvr_adj-old.rds") # Controlled by just old
if (file.exists(vols.rds)) {
  vols.dt <- readRDS(vols.rds)
} else {
  here("code/calc_hvr.R") |> source()
  #vols.dt <- vols_all.dt # Youth & Old
  vols.dt <- vols_old.dt # Just Old
}
rm(vols.rds)

## Assign VISIT labels to vols.dt
visits.lst <- list(
  VM00 = vols.dt[
    , .SD[which.min(SCANDATE)], PTID
  ][
    , .(PTID, SCANDATE, VISIT = "VM00")
  ]
)

visits.lst[["VM06"]] <- vols.dt[
  !visits.lst[["VM00"]], on = .(PTID, SCANDATE)
][
  demog.dt[VISIT == "VM06", PTID], on = "PTID", .SD[which.min(SCANDATE)], PTID
][
  , .(PTID, SCANDATE, VISIT = "VM06")
]

visits.lst[["VM12"]] <- vols.dt[
  !visits.lst[["VM00"]], on = .(PTID, SCANDATE)
][
  !visits.lst[["VM06"]], on = .(PTID, SCANDATE)
][
  demog.dt[VISIT == "VM12", PTID], on = "PTID", .SD[which.min(SCANDATE)], PTID
][
  , .(PTID, SCANDATE, VISIT = "VM12")
]

visits.lst[["VM24"]] <- vols.dt[
  !visits.lst[["VM00"]], on = .(PTID, SCANDATE)
][
  !visits.lst[["VM06"]], on = .(PTID, SCANDATE)
][
  !visits.lst[["VM12"]], on = .(PTID, SCANDATE)
][
  demog.dt[VISIT == "VM24", PTID], on = "PTID", .SD[which.min(SCANDATE)], PTID
][
  , .(PTID, SCANDATE, VISIT = "VM24")
]

visits.lst[["VM36"]] <- vols.dt[
  !visits.lst[["VM00"]], on = .(PTID, SCANDATE)
][
  !visits.lst[["VM06"]], on = .(PTID, SCANDATE)
][
  !visits.lst[["VM12"]], on = .(PTID, SCANDATE)
][
  !visits.lst[["VM24"]], on = .(PTID, SCANDATE)
][
  demog.dt[VISIT == "VM36", PTID], on = "PTID", .SD[which.min(SCANDATE)], PTID
][
  , .(PTID, SCANDATE, VISIT = "VM36")
]

visits.dt   <- rbindlist(visits.lst)

vols.dt     <- visits.dt[vols.dt, on = .(PTID, SCANDATE)]
rm(visits.dt, visits.lst)

# Missing subject is MRT62, VM24
vols.dt[is.na(VISIT), `:=`(PTID = "MRT62", VISIT = "VM24")]

# Imaging
imag.dt     <- pet.dt[vols.dt, on = .(PTID, VISIT)]

# COVARS
# Repeated MRT63 with NAs
covars.dt   <- neuropsy.dt[
  !(PTID == "MRT63" & is.na(EVALDATE))
][
  imag.dt, on = .(PTID, VISIT)
]

# All baseline data
triad.dt    <- demog.dt[covars.dt, on = .(PTID, VISIT)]
#rm(vols.dt, pet.dt, imag.dt, neuropsy.dt, covars.dt, demog.dt)

## Data cleaning
# Sex
triad.dt[
  , TAU_braak_stage := as.numeric(TAU_braak_stage)
][
  , let(
    SEX = factor(SEX, labels = c("Female", "Male")),
    # Time differences
    AGE_scan = interval(ymd(DOB), ymd(SCANDATE)) / years(1),
    EVAL_delay = ymd(EVALDATE) - ymd(SCANDATE),
    # Braak staging
    TAU_braak_group = fcase(
      TAU_braak_stage == 0, "0",
      TAU_braak_stage %in% 1:2, "1 & 2",
      TAU_braak_stage %in% 3:4, "3 & 4",
      TAU_braak_stage %in% 5:6, "5 & 6"
    )
  )
]


# MRT62 lacking SCANDATE use EVALDATE
triad.dt[
  is.na(AGE_scan),
  let(
    AGE_scan = interval(ymd(DOB), ymd(EVALDATE)) / years(1),
    EVAL_delay = 0
  )
]

# Clean DX
triad.dt[DX == "Unknown", DX := NA]
triad.dt[, DX_clean := DX]
triad.dt[DX %in% c("Atypical Dementia", "FTD", NA), DX_clean := "Other"]

# Remove useless columns
triad.dt[, c("DOB", "SCANDATE", "EVALDATE") := NULL]

# Export
saveRDS(triad.dt, here("data/rds/triad.rds"))

# TODO: Decide if remove NAs
# Must have: Imaging
#triad.dt    <- triad.dt[!is.na(AMYLOID) & !is.na(TAU_braak1) & !is.na(HVR_l)]
amy_subs.dt <- cerebra.dt[
  !is.na(AMYLOID_norm),
  .(PTID_VISIT = paste(PTID, VISIT, sep = "_"))
] |> unique()
triad.dt[, PTID_VISIT := paste(PTID, VISIT, sep = "_")]
triad.dt    <- triad.dt[amy_subs.dt, on = "PTID_VISIT"]

# Must have: Neuropsy
#triad.dt    <- triad.dt[!is.na(RAVLT_raw) & !is.na(MOCA_score)]

# No of evals
sessn       <- triad.dt[, .(SESS = as.character(.N)), PTID]
triad_bl.dt <- triad.dt[VISIT == "VM00"][sessn, on = "PTID"]
rm(sessn)

## Table1
fname       <- here("data/derivatives/table1_dx.docx")
if (!file.exists(fname) | REDOTABLES) {
  #triad_bl.dt[!is.na(APOE_n) & !is.na(MOCA_score) & DX %in% c("CN", "MCI"),
              #.(DX_clean, SEX, AGE_scan, EDUC, APOE = factor(APOE_n),
  triad_bl.dt[!is.na(MOCA_score) & !is.na(HVR_l) & DX %in% c("CN", "MCI"),
              .(DX_clean = factor(DX_clean, labels = c("NC", "MCI")),
                SEX, AGE_scan, EDUC, MOCA_score, SESS,
                #RAVLT_intro, RAVLT_raw, RAVLT_rep,
                TAU_braak_stage = TAU_braak_group,
                #AMYLOID, TAU_sum = (TAU_braak1 + TAU_braak2 + TAU_braak3 +
                                    #TAU_braak4 + TAU_braak5 + TAU_braak6),
                #HCv_l, HCv_r, HVR_l, HVR_r)] |>
                HVR = (HVR_l + HVR_r) /2)] |>
    tbl_summary(by = DX_clean,
                label = list(SEX ~ "Sex",
                             AGE_scan ~ "Age (years)",
                             EDUC ~ "Education (years)",
                             #APOE ~ "APOE4 alleles",
                             MOCA_score ~ "MoCA score",
                             #RAVLT_raw ~ "RAVLT (raw score)",
                             #RAVLT_intro ~ "RAVLT (intro score)",
                             #RAVLT_rep ~ "RAVLT (rep score)",
                             SESS ~ "Number of visits",
                             TAU_braak_stage ~ "Braak Stage (Tau)"),
                             #AMYLOID ~ "Amyloid (PET)"),
                             #TAU_sum ~ "Tau (PET)"),
                             #HCv_l ~ "HC vol (left)",
                             #HCv_r ~ "HC vol (right)",
                             #HVR_l ~ "HVR (left)",
                             #HVR_r ~ "HVR (right)"),
                statistic = all_continuous() ~ "{mean} ({sd})",
                missing_text = "Missing") |>
                #missing = "no") |>
    modify_header(label ~ "**Variable**") |>
    #modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Clinical Label**") |>
    #add_n() |>
    add_p() |> as_flex_table() |>
    flextable::save_as_docx(path = fname)
}

#fname       <- here("data/derivatives/table1_braak.docx")
#if (!file.exists(fname) | REDOTABLES) {
  #triad.dt[!is.na(APOE_n) & !is.na(MOCA_score) & !is.na(TAU_braak_group),
           #.(TAU_braak_group, SEX, AGE_scan, EDUC, APOE = factor(APOE_n),
             #MOCA_score, AMYLOID,
             ##TAU_sum = (TAU_braak1 + TAU_braak2 + TAU_braak3 +
                           ##TAU_braak4 + TAU_braak5 + TAU_braak6),
             ##HCv_l, HCv_r, HVR_l, HVR_r)] |>
             #HVR = (HVR_l + HVR_r) / 2)] |>
    #tbl_summary(by = TAU_braak_group,
                #label = list(SEX ~ "Sex",
                             #AGE_scan ~ "Age (years)",
                             #EDUC ~ "Education (years)",
                             #APOE ~ "APOE4 alleles",
                             #MOCA_score ~ "MoCA score",
                             ##RAVLT_raw ~ "RAVLT (raw score)",
                             ##RAVLT_intro ~ "RAVLT (intro score)",
                             ##RAVLT_rep ~ "RAVLT (rep score)",
                             ##SESS ~ "Number of visits",
                             ##TAU_sum ~ "Tau (PET)",
                             #AMYLOID ~ "Amyloid (PET)"),
                             ##HCv_l ~ "HC vol (left)",
                             ##HCv_r ~ "HC vol (right)",
                             ##HVR_l ~ "HVR (left)",
                             ##HVR_r ~ "HVR (right)"),
                #statistic = all_continuous() ~ "{mean} ({sd})",
                ##missing_text = "Missing") |>
                #missing = "no") |>
    #modify_header(label ~ "**Variable**") |>
    ##modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Clinical Label**") |>
    ##add_n() |>
    #add_p() |> as_flex_table() |>
    #flextable::save_as_docx(path = fname)
#}
