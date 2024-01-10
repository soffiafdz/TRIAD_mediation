#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(lubridate)
library(stringr)

## Dependencies

# Load demographic data for subsetting
if (!exists("demog.dt")) {
  here("code/demographics.R") |> source()
}

hcvc.dt     <- here("data/derivatives/vols_hcvc.csv") |> fread()
hcvc.dt[, `:=`(PTID = str_extract(ID, "(?<=stx_).*(?=_20)"),
               SCANDATE = ymd(str_extract(ID, "\\d{4}-\\d{2}-\\d{2}")))]
hcvc.dt[is.na(PTID), PTID := "MRT62"]

scale.dt    <- here("data/derivatives/icc_scale.csv") |> fread()
scale.dt[, SCANDATE := ymd(SCANDATE)]

# Parse HCvols on native space
vols.dt     <- hcvc.dt[scale.dt, on = .(PTID, SCANDATE),
                       .(PTID, SCANDATE,
                         ICC    = ICC  / 1000,
                         HC_l   = LHC  / (SCALEFACTOR * 1000),
                         HC_r   = RHC  / (SCALEFACTOR * 1000),
                         CSF_l  = LCSF / (SCALEFACTOR * 1000),
                         CSF_r  = RCSF / (SCALEFACTOR * 1000))]
rm(hcvc.dt, scale.dt)

vols.dt     <- melt(vols.dt, id.vars = c("PTID", "SCANDATE", "ICC"),
                    variable.name = "ROI", value.name = "VAL")

# Adjust for ICV
# Linear model of VAL ~ ROI & average ICC on healthy people
# Young & Old
cn1.dt      <- vols.dt[PTID %in% demog.dt[DX %in% c("CN", "Young"),
                                          unique(PTID)],
                       .SD[which.min(SCANDATE)],
                       .(PTID, ROI)]
all_icc     <- cn1.dt[!duplicated(PTID), mean(ICC)]
b1.dt       <- cn1.dt[, .(b_all = summary(lm(VAL ~ ICC))$coefficients[2]), ROI]

# Just Old
cn2.dt      <- vols.dt[PTID %in% demog.dt[DX == "CN", unique(PTID)],
                       .SD[which.min(SCANDATE)],
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
vols_all.dt <- vols.dt[, .(PTID, SCANDATE, ROI, VAL_adj_all)] |>
  dcast(... ~ ROI)
vols_all.dt <- vols_all.dt[, .(PTID, SCANDATE, HCv_l = HC_l, HCv_r = HC_r,
                               HVR_l = HC_l / (HC_l + CSF_l),
                               HVR_r = HC_r / (HC_r + CSF_r))]
write_rds(vols_all.dt, here("data/rds/hcv_hvr_adj-all.rds"))

vols_old.dt <- vols.dt[, .(PTID, SCANDATE, ROI, VAL_adj_old)] |>
  dcast(... ~ ROI)
vols_old.dt <- vols_old.dt[, .(PTID, SCANDATE, HCv_l = HC_l, HCv_r = HC_r,
                               HVR_l = HC_l / (HC_l + CSF_l),
                               HVR_r = HC_r / (HC_r + CSF_r))]
write_rds(vols_old.dt, here("data/rds/hcv_hvr_adj-old.rds"))
rm(vols.dt)
