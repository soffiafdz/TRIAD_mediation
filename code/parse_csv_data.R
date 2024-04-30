#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
library(stringr)
library(lubridate)

## Read files
# File names
fname1  <- here("data/TRIAD_data2024.csv")
fname2  <- here("data/TRIAD_Raket_all.csv")
fname3  <- here("data/derivatives/wmh_vols.csv")

# Parsing
DT1     <- fread(fname1)
DT2     <- fread(fname2)
DT3     <- fread(fname3)
rm(fname1, fname2, fname3)

## DTs
# Demographics
cols    <- c("id", "visit", "dx", "dob", "gender", "education", "apoe",
             "cdr", "mmse",
             "date_mri", "mri_qc", "date_mk", "mk_qc", "date_nav", "nav_qc")
DT1     <- DT1[, ..cols]

setnames(DT1, cols,
         c("PTID", "VISIT", "DX", "DOB", "SEX", "EDUC", "APOE",
           "CDR", "MMSE",
           "DATE_mri", "QC_mri", "DATE_mk", "QC_mk", "DATE_nav", "QC_nav"))

rm(cols)

excl.dt  <- DT1[DT1[, .N, .(PTID, DATE_mri)][N != 1],
               on = .(PTID, DATE_mri)
               ][DATE_mk == "" | DATE_nav == "",
               .(PTID, VISIT)]

DT1     <- DT1[!excl.dt, on = .(PTID, VISIT)]

## Keep only certain DX categories
dx_keep <- c("CN", "CN(Y)", "CN (Y)",
             "SCI", "SCI (mixed vascular component)",
             "MCI", "MCI not due to AD", "early MCI", "MCI (mixed)", "aMCI",
             "AD", "early AD", "early onset AD",
             "mild AD", "Mild AD",
             "EOAD", "Possible EOAD",
             "AD (with an atypical distribution of Amyloid-beta)")

#DT1     <- DT1[DX %in% dx_keep]

# Manage dates
DT1[, `:=`(DOB      = dmy(DOB),
           DATE_mri = dmy(DATE_mri),
           DATE_mk  = dmy(DATE_mk),
           DATE_nav = dmy(DATE_nav))]

# Calculate age (using dob and midpoint between earliest and latest scan)
DT1[, `:=`(DATE_earliest  = do.call(pmin, .SD),
           DATE_latest    = do.call(pmax, .SD)),
    .SDcols = grep("DATE", names(DT1), value = TRUE)]

DT1[, DATE_midpoint := DATE_earliest + (DATE_latest - DATE_earliest) / 2]

DT1[, AGE := as.period(interval(DOB, DATE_midpoint))$year]

# Sex
DT1[, SEX := factor(SEX, labels = c("Female", "Male"))]

## Remove people younger than 40 (?)
covars.dt     <- covars.dt[AGE >= 40]

write_rds(DT1, here("data/rds/covars.rds"))


# Raket disease offset
cols    <- c("RID", "visit", "nav_adni_suvr_fullcg_neocortex", "AB", "EDT")
DT2     <- DT2[, ..cols]

setnames(DT2, cols,
         c("PTID", "VISIT", "AB_neocortex", "AB_bool", "RAKET_edt"))

rm(cols)

## AB positivity
## nav_adni_suvr_fullcg_neocortex > 1.55
# Separate groups by Raket disease offset
#DT2[RAKET_edt < -1,                    RAKET_group := "<-1"  ]
#DT2[RAKET_edt >= -1  & RAKET_edt < 0,  RAKET_group := "-1—0" ]
#DT2[RAKET_edt >= 0   & RAKET_edt < 2,  RAKET_group := "0—2"  ]
#DT2[RAKET_edt >= 2   & RAKET_edt < 4,  RAKET_group := "2—4"  ]
#DT2[RAKET_edt >= 4,                    RAKET_group := ">4"   ]

DT2[RAKET_edt <= 0, RAKET_group := 0]
DT2[RAKET_edt > 0 & RAKET_edt <= 3, RAKET_group := 1]
DT2[RAKET_edt > 3, RAKET_group := 2]

DT2[, RAKET_group := factor(RAKET_group,
                            levels = 0:2,
                            labels = c("Healthy",
                                       "Early stages",
                                       "Late stages"))]

write_rds(DT2, here("data/rds/raket_eds.rds"))

# WMH
setnames(DT3, "SESSION", "VISIT")
write_rds(DT3, here("data/rds/wmh_vols.rds"))

rm(DT1, DT2, DT3)
