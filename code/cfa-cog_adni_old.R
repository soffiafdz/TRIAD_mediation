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

cols  <- c("PTID", "VISCODE", "EXAMDATE", "DX", COI)
DT    <- adnimerge[, ..cols]

## Abbreviations
#COI_a <- c("A11", "A13", "AQ4", "CDRSB", "MMSE", "PACC1", "PACC2",
           #"REYf", "REYi", "REYl", "REYp")

COI_used    <- c("ADASQ4", "CDRSB", "MMSE")
COI_used2   <- c("ADASQ4", "PACC.digits", "RAVLT.learning")

setnames(DT, COI[6], COI_used2[2])

### EFA
## Large sample (COI1): Best estimates (CFI > .9) are for 4,5 factors
#nf    <- 1:5
#fpath <- here("data/rds/adni_efa-cog1.rds")
#if (!file.exists(fpath) | REFIT_FA) {
  #cols      <- c("DX", COI_used)
  #Dxs       <- DT[, levels(DX)]
  #efa.f     <- vector("list", length(Dxs))
  #for (i in seq_along(Dxs)) {
    #efa.f[[i]] <- DT[Dxs[i], on = "DX", ..COI_used] |>
    #na.omit() |>
    #efa(nfactors = nf)
  #}
  #saveRDS(efa.f)
#} else {
  #efa.f     <- readRDS(fpath)
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
CogDecline =~ ADASQ4 + CDRSB + MMSE
!f2 =~ REYi + REYl
'

cfa2.m      <- '
CogDecline =~ ADASQ4 + PACC.digits + RAVLT.learning
!f2 =~ REYi + REYl
'

#fpath <- here("data/rds/adni_cfa-cog1.rds")
#if (!file.exists(fpath) | REFIT_FA) {
  #cols      <- c("DX", COI_used)
  #cfa1.f    <- DT[, ..cols] |>
  #na.omit() |>
  #cfa(model = cfa.m, group = "DX", meanstructure = T) |>
  #saveRDS(fpath)
#} else {
  #cfa1.f    <- readRDS(fpath)
#}

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
  #saveRDS(fpath)
#} else {
  #cfa2.f     <- readRDS(fpath)
#}

## Metric invariance: Only with F1
#fpath <- here("data/rds/adni_cfa-cog2.rds")
#if (!file.exists(fpath) | REFIT_FA) {
  #cols      <- c("DX", COI_used)
  #cfa2.f    <- DT[, ..cols] |>
  #na.omit() |>
  #cfa(model = cfa.m, group = "DX", meanstructure = T,
      #group.equal = c("loadings")) |>
  #saveRDS(fpath)
#} else {
  #cfa2.f     <- readRDS(fpath)
#}

### Not comparing between groups, though.
### Some CNs will turn into AD.

#fpath <- here("data/rds/adni_cfa-cog.rds")
#if (!file.exists(fpath) | REFIT_FA) {
  ##cols      <- c(COI_used)
  #cfa.f     <- DT[, ..COI_used] |>
  #na.omit() |>
  #cfa(model = cfa.m, meanstructure = T) |>
  #saveRDS(fpath)
#} else {
  #cfa.f     <- readRDS(fpath)
#}

#fpath <- here("data/rds/adni_cfa-cog2.rds")
#if (!file.exists(fpath) | REFIT_FA) {
  #cols      <- c("PTID", COI_used2)
  #cfa2.f    <- DT[, ..cols] |>
  #na.omit() |>
  #cfa(model = cfa2.m, meanstructure = T, cluster = "PTID", std.lv = T) |>
  #saveRDS(fpath)
#} else {
  #cfa2.f     <- readRDS(fpath)
#}

## Do final CFA on the data with MRI
cols        <- c("PTID", "VISCODE", "EXAMDATE", COI_used, COI_used2) |> unique()
DT_mri      <- DT[dx.dt, on = .(PTID, VISCODE)][, ..cols] |> na.omit()
fpath <- here("data/rds/adni_cfa-cog_mri.rds")
if (!file.exists(fpath) | REFIT_FA) {
  #cols      <- c(COI_used)
  cfa_m.f   <- cfa(cfa.m, DT_mri, meanstructure = T, std.lv = T)
  saveRDS(cfa_m.f, fpath)
} else {
  cfa_m.f   <- readRDS(fpath)
}

## Fit MCFA: cluster on SUB
fpath <- here("data/rds/adni_mcfa-cog_mri.rds")
if (!file.exists(fpath) | REFIT_FA) {
  #cols      <- c(COI_used)
  mcfa_m.f  <- cfa(cfa.m, DT_mri, meanstructure = T, std.lv = T,
                   cluster = "PTID")
  saveRDS(mcfa_m.f, fpath)
} else {
  mcfa_m.f   <- readRDS(fpath)
}

## Fit MCFA: cluster on SUB
fpath <- here("data/rds/adni_mcfa-cog_mri2.rds")
if (!file.exists(fpath) | REFIT_FA) {
  #cols      <- c(COI_used)
  mcfa_m2.f  <- cfa(cfa2.m, DT_mri, meanstructure = T, std.lv = T,
                   cluster = "PTID")
  saveRDS(mcfa_m2.f, fpath)
} else {
  mcfa_m2.f   <- readRDS(fpath)
}

outmds       <- here('data/derivatives',
                     paste0('mcfa-cognition_mri', 1:2, '.md'))

for (outmd in outmds) {
  sink(outmd)
  parameters::parameters(mcfa_m.f) |> insight::print_md()
  sink()
  rm(outmd)
}

## Save Cognition's latent factor
outrds      <- here('data/rds/adni_cog-latent.rds')
DT_mri[, COG_latent := lavPredict(mcfa_m.f)]
DT_mri[, COG_latent2 := lavPredict(mcfa_m2.f)]
cog.dt      <- DT_mri
saveRDS(cog.dt, outrds)
