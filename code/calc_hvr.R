#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(lubridate)
library(stringr)

## Dependencies
hcvc    <- here("data/derivatives/vols_hcvc.csv") |> fread()
hcvc[, `:=`(PTID = str_extract(ID, "(?<=stx_).*(?=_20)"),
                 SCANDATE = ymd(str_extract(ID, "\\d{4}-\\d{2}-\\d{2}")))]

vols.dt <- hcvc[, .(PTID, SCANDATE,
                    HCv_l = LHC / 1000,
                    HCv_r = RHC / 1000,
                    HVR_l = LHC / (LHC + LCSF),
                    HVR_r = RHC / (RHC + RCSF))]

write_rds(vols.dt, here("data/rds/hcv_hvr.rds"))
rm(hcvc)
