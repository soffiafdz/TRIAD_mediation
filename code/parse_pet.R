#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)

## Read/Parse CSV files
fpaths      <- here("data/rds",
                    c("native_nav_cerebra.rds", "native_mk_cerebra.rds"))

if (any(!file.exists(fpaths))) here("code/parse_csv_data.R") |> source()

nav.dt      <- read_rds(fpaths[1]) |> setkey(PTID, VISIT, LABEL_id)
mk.dt       <- read_rds(fpaths[2]) |> setkey(PTID, VISIT, LABEL_id)

rm(fpaths)

## Amyloid & Tau by ROI (CEREBRA)
# Native
nav_lst <- here("data/data_2024") |>
  list.files(pattern = "native_nav", recursive = TRUE, full.names = TRUE)

nav.dt  <- nav_lst |>
  lapply(fread) |>
  rbindlist(idcol = TRUE)

nav_lab <- nav_lst |>
  str_extract("(?<=native_nav_).*(?=_cerebra.csv)") |>
  data.table()

nav_lab[, `:=`(.id    = 1:.N,
               PTID   = str_split_i(V1, "_", 1),
               VISIT  = str_split_i(V1, "_", 2))]

nav.dt  <- nav_lab[, -"V1"
                   ][nav.dt[, -c("mx", "my", "mz")], on = ".id"
                   ][, -".id"]

rm(nav_lst, nav_lab)

# Normalize SUVR values by Avg L/R cerebellar gray matter
nav.dt     <- nav.dt[id != "label"]
nav.dt     <- nav.dt[id %like% "Cerebellum_Gray_Matter",
               .(CGM = mean(val)),
               .(PTID, VISIT)
               ][nav.dt, on = .(PTID, VISIT),
               .(PTID, VISIT, LABEL_id = id, VOLUME_nav = volume,
                 SUVR_nav = val / CGM)]

# Normalize SUVR values by ROIvol
nav.dt[, SUVR_nav_norm := SUVR_nav / VOLUME_nav]

write_rds(nav.dt, here("data/rds/native_nav_cerebra.rds"))

# MK
mk_lst  <- here("data/data_2024") |>
  list.files(pattern = "native_mk", recursive = TRUE, full.names = TRUE)

mk.dt     <- mk_lst |>
  lapply(fread) |>
  rbindlist(idcol = TRUE)

mk_lab  <- mk_lst |>
  str_extract("(?<=native_mk_).*(?=_cerebra.csv)") |>
  data.table()

mk_lab[, `:=`(.id    = 1:.N,
              PTID   = str_split_i(V1, "_", 1),
              VISIT  = str_split_i(V1, "_", 2))]

mk.dt     <- mk_lab[, -"V1"
                  ][mk.dt[, -c("mx", "my", "mz")], on = ".id"
                  ][, -".id"]

rm(mk_lst, mk_lab)

# Normalize SUVR values by Avg L/R cerebellar gray matter
mk.dt     <- mk.dt[id != "label"]
mk.dt     <- mk.dt[id %like% "Cerebellum_Gray_Matter",
               .(CGM = mean(val)),
               .(PTID, VISIT)
               ][mk.dt, on = .(PTID, VISIT),
               .(PTID, VISIT, LABEL_id = id, VOLUME_mk = volume,
                 SUVR_mk = val / CGM)]

# Normalize SUVR values by ROIvol
mk.dt[, SUVR_mk_norm := SUVR_mk / VOLUME_mk]

write_rds(mk.dt, here("data/rds/native_mk_cerebra.rds"))

# Keep only Visits with both Amyloid and Tau
pet.dt      <- mk.dt[nav.dt, on = .(PTID, VISIT, LABEL_id)]
rm(mk.dt, nav.dt)

## Export object
write_rds(pet.dt, here("data/rds/pet_cerebra.rds"))
