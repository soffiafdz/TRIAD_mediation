#!/usr/bin/env Rscript

library(here)
library(data.table)

COI         <- c("SubjectID", "VisitID", "WMH_vol")

corr_date.f <- function(x) x[, VisitID := lubridate::ymd(VisitID)]

DT          <- here("data/adni/vols") |>
list.files(full.names = TRUE) |>
lapply(fread, select = COI) |>
lapply(corr_date.f) |>
rbindlist() |>
setnames(c("PTID", "EXAMDATE", "WMHvol")) |>
fwrite(here("data/adni_wmh-vol.csv"))
