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

## IN
fpath       <- here("data/adni_raket_input_neda.csv")
if (file.exists(fpath)) {
  DT        <- fread(fpath)
} else {
  sprintf("File: %s is required but could not be found.", fpath) |> stop()
}
## Deprecated ##
## Keep only A+
#cog.dt      <- DT[.(1), on = "AmyloidStatus", c(1, 4, 18:29)] |>

## Concatenate AB- CN & AB+ CN+MCI+AD
cog.dt      <- DT[DX_new == "CNNEG" | AmyloidStatus,
                  .SD,
                  by = .(RID, EXAMDATE),
                  .SDcols = MMSE:MOCA] |>
na.omit()

## CFA with 4 factors
## This didn't work.
## Left here as evidenced of failed work.
#cfa.m       <- '
#COG =~ MOCA + MMSE + ADAS13 + mPACCtrailsB
#MEM =~ ADASQ4 + RAVLT.learning + RAVLT.forgetting + mPACCdigit
#GEN =~ COG + MEM
#'

models.lst  <- list(COG1 = 'COG =~ MOCA + ADAS13 + MMSE',
                    COG2 = 'COG =~ MOCA + ADAS13 + CDRSB')

## Fit MCFA: cluster on SUB
mcfa.lst <- cogl.lst <- vector("list", length(models.lst))
paths       <- list(rds = "data/rds", der = "data/derivatives")
for (i in seq_along(models.lst)) {
  fnames    <- list(mod = sprintf("adni_raket_mcfa-cog_mod%i", i),
                    var = sprintf("adni_raket_mcfa-cog_var%i", i))

  # Model fitting
  fpath     <- sprintf("%s/%s.%s", paths$rds, fnames$mod, "rds") |> here()
  if (all(!REFIT_FA, file.exists(fpath))) {
    mcfa.lst[[i]] <- readRDS(fpath)
  } else {
    mcfa.lst[[i]] <- cfa(models.lst[[i]], cog.dt,
                         meanstructure  = TRUE,     # To compare models
                         std.lv         = TRUE,     # Standardized coefficients
                         cluster        = "RID")    # Random effects for Subj
    saveRDS(mcfa.lst[[i]], fpath)
  }

  ### This has to be done in the terminal, doesn't work in Rscript
  ## Printing parameters
  #fpath     <- sprintf("%s/%s.%s", paths$der, fnames$mod, "md") |> here()
  #if (any(REFIT_FA, !file.exists(fpath))) {
    #sink(fpath)
    #parameters::parameters(mcfa.lst[[i]]) |> insight::print_md()
    #sink()
  #}

  # Obtain latent cognitive measures
  cogl.lst[[i]] <- cog.dt[, .(RID, EXAMDATE, lavPredict(mcfa.lst[[i]]))]

  ## Save RDS of latent variable
  fpath     <- sprintf("%s/%s.%s", paths$rds, fnames$var, "rds") |> here()
  if (any(REFIT_FA, !file.exists(fpath))) {
    cogl.lst[[i]][cog.dt, on = .(RID, EXAMDATE)] |> saveRDS(fpath)
  }

  # Save CSV of latent variable
  fpath     <- sprintf("%s/%s.%s", paths$der, fnames$var, "csv") |> here()
  if (any(REFIT_FA, !file.exists(fpath))) {
    cogl.lst[[i]][cog.dt, on = .(RID, EXAMDATE)] |> fwrite(fpath)
  }
  rm(fpath)
}
