#!/usr/bin/env Rscript

library(here)
library(data.table)
library(readr)
library(ADNIMERGE)
library(lavaan)

data(adnimerge)
setDT(adnimerge)

## Imputed Dxs from ADNIMERGE
REFIT_FA    <- T
USE_IMPUTED <- T

fpath       <- here("data/rds/adni_dxs_imputed.rds")
if (file.exists(fpath)) {
  dx.dt     <- readRDS(fpath)
} else {
  here("code/impute_dx.R") |> source()
}

if (USE_IMPUTED) dx.dt[, DX := stringr::str_remove(DX, "\\?")]

# Covariates of interest
# N ~11k & ~ 6k on MRIdataset
COI1  <- c("ADAS", "CDRSB", "MMSE", "PACC", "RAVLT")
# N ~6k & ~2.5K on MRIdataset
COI2  <- c("Ecog","MOCA")

## Will focus only on COI1. Leave this if needed
#COI   <- c(COI1, COI2) |>
#paste(collapse = "|") |>
#grep(x = names(adnimerge), value = TRUE) |>
#grep(pattern = "bl", invert = TRUE, value = TRUE) |>
#sort()

#COI1  <- COI1 |> paste(collapse = "|") |> grep(x = COI, value = TRUE)
#COI2  <- COI2 |> paste(collapse = "|") |> grep(x = COI, value = TRUE)


COI   <- COI1 |>
paste(collapse = "|") |>
grep(x = names(adnimerge), value = TRUE) |>
grep(pattern = "bl", invert = TRUE, value = TRUE) |>
sort()

cols  <- c("PTID", "VISCODE", "DX", COI)
DT    <- adnimerge[, ..cols]

## Abbreviations
COI_a <- c("A11", "A13", "AQ4", "CDRSB", "MMSE", "PACC1", "PACC2",
           "REYf", "REYi", "REYl", "REYp")

setnames(DT, COI, COI_a)

COI_used    <- c("AQ4", "CDRSB", "MMSE")

### EFA
## Large sample (COI1): Best estimates (CFI > .9) are for 4,5 factors
#nf    <- 1:5
#fpath <- here("data/rds/efa-cog1_adni.rds")
#if (!file.exists(fpath) | REFIT_FA) {
  #cols      <- c("DX", COI_used)
  #Dxs       <- DT[, levels(DX)]
  #efa.f     <- vector("list", length(Dxs))
  #for (i in seq_along(Dxs)) {
    #efa.f[[i]] <- DT[Dxs[i], on = "DX", ..COI_used] |>
    #na.omit() |>
    #efa(nfactors = nf)
  #}
  #write_rds(efa.f)
#} else {
  #efa.f     <- read_rds(fpath)
#}
##summary(efa.f$nf4)

# CFA with 4 factors
#cfa.m       <- '
#f1 =~ A11 + A13 + CDRSB + REYi
#f2 =~ A13 + AQ4 + MMSE + PACC1 + PACC2 + REYf + REYi
#f3 =~ A13 + CDRSB + MMSE + PACC1 + REYf + REYl + REYp
#f4 =~ MMSE + PACC1 + REYf + REYl + REYp
#f5 =~ A13 + MMSE + PACC1 + PACC2 + REYi + REYl + REYp
#'

## 2Factor model reaches configural invariance
## Removing the memory factor reaches metric invariance
cfa.m       <- '
f1 =~ AQ4 + CDRSB + MMSE
!f2 =~ REYi + REYl
'
fpath <- here("data/rds/cfa-cog1_adni.rds")
if (!file.exists(fpath) | REFIT_FA) {
  cols      <- c("DX", COI_used)
  cfa1.f    <- DT[, ..cols] |>
  na.omit() |>
  cfa(model = cfa.m, group = "DX", meanstructure = T) |>
  write_rds(fpath)
} else {
  cfa1.f    <- read_rds(fpath)
}

## Second order factor has negative variance.
## Keep one factor for cognition and another for memory
#cfa2.m      <- '
#f1 =~ A13 + CDRSB + MMSE
#f2 =~ REYi + REYl
#f3 =~ 1*f1 + 1*f2
#f3 ~~ f3
#'
#if (!file.exists(fpath) | REFIT_FA) {
  #cols      <- c("DX", COI_used)
  #cfa2.f    <- DT[, ..cols] |>
  #na.omit() |>
  #cfa(model = cfa2.m, group = "DX") |>
  #write_rds(fpath)
#} else {
  #cfa2.f     <- read_rds(fpath)
#}

## Metric invariance: Only with F1
fpath <- here("data/rds/cfa-cog2_adni.rds")
if (!file.exists(fpath) | REFIT_FA) {
  cols      <- c("DX", COI_used)
  cfa2.f    <- DT[, ..cols] |>
  na.omit() |>
  cfa(model = cfa.m, group = "DX", meanstructure = T,
      group.equal = c("loadings")) |>
  write_rds(fpath)
} else {
  cfa2.f     <- read_rds(fpath)
}

### Not comparing between groups, though.
### Some CNs will turn into AD.

fpath <- here("data/rds/cfa-cog_adni.rds")
if (!file.exists(fpath) | REFIT_FA) {
  #cols      <- c(COI_used)
  cfa.f     <- DT[, ..COI_used] |>
  na.omit() |>
  cfa(model = cfa.m, meanstructure = T) |>
  write_rds(fpath)
} else {
  cfa.f     <- read_rds(fpath)
}


fpath <- here("data/rds/cfa-cog_mri_adni.rds")
if (!file.exists(fpath) | REFIT_FA) {
  #cols      <- c(COI_used)
  cfa_m.f    <- DT[dx.dt, on = .(PTID, VISCODE)][, ..COI_used] |>
  na.omit() |>
  cfa(model = cfa.m, meanstructure = T, std.lv = T) |>
  write_rds(fpath)
} else {
  cfa_m.f    <- read_rds(fpath)
}
