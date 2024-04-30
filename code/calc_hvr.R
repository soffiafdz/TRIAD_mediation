#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(lubridate)
library(stringr)

## Dependencies
fpaths      <- here("data/derivatives", c("icc_scale.csv", "vols_hcvc.csv"))
scale.dt    <- fpaths[1] |> fread()
hcvc.dt     <- fpaths[2] |> fread()

hcvc.dt[, `:=`(PTID   = str_extract(ID, "(?<=stx_).*(?=_VM)"),
               VISIT  = str_extract(ID, "VM\\d*(?=_t1)"))]

# Load demographic data for subsetting
fpaths      <- here("data/rds", c("covars.rds", "raket_eds.rds"))
if (any(!file.exists(fpaths)))  here("code/parse_csv_data.R") |> source()

covars.dt   <- read_rds(fpaths[1]) |> setkey(PTID, VISIT)
raket.dt    <- read_rds(fpaths[2]) |> setkey(PTID, VISIT)

subset.dt   <- covars.dt[, .(PTID, VISIT, DX, DATE_mri)]

subset.dt   <- subset.dt[raket.dt[AB_bool == TRUE, .(PTID, VISIT)]]

rm(covars.dt, raket.dt)

# Parse HCvols on native space
vols.dt     <- scale.dt[hcvc.dt, on = .(PTID, VISIT),
                       .(PTID, VISIT,
                         ICC    = ICC  / 1000,
                         HC_l   = LHC  / (SCALEFACTOR * 1000),
                         HC_r   = RHC  / (SCALEFACTOR * 1000),
                         CSF_l  = LCSF / (SCALEFACTOR * 1000),
                         CSF_r  = RCSF / (SCALEFACTOR * 1000))]
rm(hcvc.dt, scale.dt)

vols.dt     <- melt(vols.dt, id.vars = c("PTID", "VISIT", "ICC"),
                    variable.name = "ROI", value.name = "VAL")

setkey(vols.dt, PTID, VISIT)

# Adjust for ICV
# Linear model of VAL ~ ROI & average ICC on healthy people
# Young & Old
cn1.dt      <- vols.dt[subset.dt[DX == "CN" | DX %like% "Y"]
                       ][!is.na(VAL),
                       .SD[which.min(DATE_mri)],
                       .(PTID, ROI)]
all_icc     <- cn1.dt[!duplicated(PTID), mean(ICC)]
b1.dt       <- cn1.dt[, .(b_all = summary(lm(VAL ~ ICC))$coefficients[2]), ROI]

# Just Old
cn2.dt      <- vols.dt[subset.dt[DX == "CN"]
                       ][!is.na(VAL),
                       .SD[which.min(DATE_mri)],
                       .(PTID, ROI)]
old_icc     <- cn2.dt[!duplicated(PTID), mean(ICC)]
b2.dt       <- cn2.dt[, .(b_old = summary(lm(VAL ~ ICC))$coefficients[2]), ROI]

# Merge slopes
b.dt        <- b1.dt[b2.dt, on = "ROI"]
rm(cn1.dt, cn2.dt, b1.dt, b2.dt)

#vols.dt <- dcast(vols.dt, ... ~ ROI, value.var = "VAL")

# Adjust by head size
vols.dt[, `:=`(VAL_adj_all = VAL - b.dt[ROI == ROI, b_all] * (ICC - all_icc),
               VAL_adj_old = VAL - b.dt[ROI == ROI, b_old] * (ICC - old_icc))]
rm(b.dt, all_icc, old_icc)

# All
vols_all.dt <- vols.dt[, .(PTID, VISIT, ROI, VAL_adj_all)] |>
  dcast(... ~ ROI)
vols_all.dt <- vols_all.dt[, .(PTID, VISIT, HCv_l = HC_l, HCv_r = HC_r,
                               HVR_l = HC_l / (HC_l + CSF_l),
                               HVR_r = HC_r / (HC_r + CSF_r))]
write_rds(vols_all.dt, here("data/rds/hcv_hvr_adj-young-old.rds"))

vols_old.dt <- vols.dt[, .(PTID, VISIT, ROI, VAL_adj_old)] |>
  dcast(... ~ ROI)
vols_old.dt <- vols_old.dt[, .(PTID, VISIT, HCv_l = HC_l, HCv_r = HC_r,
                               HVR_l = HC_l / (HC_l + CSF_l),
                               HVR_r = HC_r / (HC_r + CSF_r))]
write_rds(vols_old.dt, here("data/rds/hcv_hvr_adj-old.rds"))
rm(vols.dt)
