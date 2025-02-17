#!/usr/bin/env Rscript

library(here)
library(data.table)
#library(readr)
#library(ADNIMERGE)
library(lavaan)

#data(adnimerge)
#setDT(adnimerge)

### CONSTANTS
REFIT_FA    <- F
USE_IMPUTED <- T

### INPUT
## Cognitive tests
tnames      <- c("NEUROBAT", "ADAS_ADNIGO23", "MMSE", "MOCA")
fpaths      <- tnames |>
  sprintf(fmt = "data/adni/cognitive_tests/%s_Jan2025.csv") |>
  here()

# Clean ADAS name
tnames      <- gsub("_.*", "", tnames)

tests.lst   <- fpaths |>
  lapply(function(fpath){
    fname   <- basename(fpath)
    if (!file.exists(fpath)) sprintf(
      "File: %s is required but could not be found",
      fpath
    ) |> stop()
    fread(fpath, key = c("PTID", "VISCODE2"))
})


names(tests.lst) <- tnames

## Imputed Dxs from ADNIMERGE
fpath       <- here("data/rds/adni_dxs_imputed.rds")
if (file.exists(fpath)) {
  dx.dt     <- readRDS(fpath)
} else {
  here("code/impute_dx.R") |> source()
}

if (USE_IMPUTED) dx.dt[, DX := stringr::str_remove(DX, "\\?")]
setkey(dx.dt, PTID, VISCODE) # In case of using visits of interest
rm(fpaths, fpath)

### Data CLEANING
## Relevant items
## See '{proj}/data/cogdomains.md'
items.lst   <- list(
  Memory = list(
    NEUROBAT = c(
      # WMS-R: Logical Memory - immediate & recall
      "LIMMTOTAL", "LDELTOTAL",
      # Rey-AVLT:  Trials 1:6
      "AVTOT1", "AVTOT2", "AVTOT3", "AVTOT4", "AVTOT5", "AVTOT6",
      # Rey-AVLT:  Trial B, 30 min delay & recognition
      "AVTOTB", "AVDEL30MIN", "AVDELTOT"
    ),
    ADAS = c(
      "Q1SCORE", # Word recall
      "Q4SCORE", # Delayed word recall
      "Q7SCORE", # Orientation
      "Q8SCORE"  # Word recognition
    ),
    MMSE = c(
      # Orientation
      "MMDATE", "MMYEAR", "MMMONTH", "MMDAY", "MMSEASON",
      "MMHOSPIT", "MMFLOOR", "MMCITY", "MMAREA", "MMSTATE",
      # Word (Ball): immediate/delayed recall
      "WORD1", "WORD1DL"
    ),
    MOCA = c(
      ## These need to be cleaned ##
      # Sum of: Trials registration
      "IMMT1W1", "IMMT1W2", "IMMT1W3", "IMMT1W4", "IMMT1W5",
      "IMMT2W1", "IMMT2W2", "IMMT2W3", "IMMT2W4", "IMMT2W5",
      # Sum of: Delayed recall of word list
      # These need to be cleaned: keep only `without cue` (1)
      "DELW1", "DELW2", "DELW3", "DELW4", "DELW5"
      ## These need to be cleaned ##
    )
  ),
  ExecFun = list(
    NEUROBAT = c(
      # Clock
      "CLOCKCIRC", "CLOCKSYM", "CLOCKNUM", "CLOCKHAND", "CLOCKTIME",
      # WAIS-R: Digit span backward
      #"DSPANBAC",
      ## WMS-R: Digit symbol & Digit span forward
      #"DIGITSCOR", "DSPANFOR"
      # Trails A/B Time
      "TRAASCOR", "TRABSCOR"
    ),
    ADAS = "Q13SCORE",    # Number cancellation task
    #MMSE = "WORLDSCORE",  # Spell `world` backwards
    MOCA = c(
      # Abstraction
      "ABSMEAS", "ABSTRAN",
      # Trails
      "TRAILS",
      # Digits backward/forward
      "DIGBACK", "DIGFOR",
      ## This needs to be cleaned ##
      # Serial 7 (Total)
      "SERIAL1", "SERIAL2", "SERIAL3", "SERIAL4", "SERIAL5",
      ## This needs to be cleaned ##
      # Letters/tapping (Errors)
      "LETTERS"
    )
  ),
  Language = list(
    NEUROBAT = c(
      # Category fluency: animals & vegetables
      "CATANIMSC"#, "CATVEGESC",
      # Boston naming test
      #"BNTTOTAL"
    ),
    ADAS = c(
      "Q2SCORE",  # Commands
      "Q5SCORE",  # Naming
      "Q6SCORE"   # Ideational-Praxis
    ),
    MMSE = c(
    "MMREPEAT",   # Repeat after instructor
    "MMHAND",     # Take paper
    "MMFOLD",     # Fold paper
    "MMONFLR",    # Place paper on floor
    "MMREAD",     # Read paper
    "MMWRITE"     # Write a sentence
    ),
    MOCA = c(
      # Naming
      "CAMEL", "LION", "RHINO",
      # Repeat sentence
      "REPEAT1", "REPEAT2",
      # Fluency-F
      "FFLUENCY"
    )
  ),
  VisospatialFun = list(
    NEUROBAT = c(
      # Clock
      "COPYCIRC", "COPYSYM", "COPYNUM", "COPYTIME"
    ),
    ADAS = "Q3SCORE", # Constructional Praxis
    MMSE = "MMDRAW"   # Interlocking pentagons
  )
)

# Create a list for all cognitive domains
cog_data.lst <- list()
for (domain in names(items.lst)) {
  # Cognitive domain
  cog_data.lst[[domain]] <- list()
  for (test in tnames) {
    cols <- c("PTID", "VISCODE2", "VISDATE", items.lst[[domain]][[test]])
    if (length(cols) > 3) {
      cog_data.lst[[domain]][[test]] <-
        ## All available data for specific subjects
        #tests.lst[[test]][
          #!is.na(VISDATE)
        #][
          #PTID %in% unique(dx.dt$PTID),
          #..cols
        #] |>
        ## Only visits of interest
         tests.lst[[test]][
           dx.dt,
           nomatch = NULL,
           ..cols
         ] |>
        melt(id = 1:3) |>
        suppressWarnings() |>
        na.omit() |>
        unique()
    }
  }
  cog_data.lst[[domain]] <- cog_data.lst[[domain]] |>
    rbindlist() |>
    dcast(... ~ variable) #|>
    ## Keep only most recent visit
    #{function(DT) DT[order(-VISDATE)][!duplicated(PTID)]}()
}
rm(cols, domain, test, tests.lst, tnames)

## Further cleaning
# Memory — MoCA: Sum of trials of registration
cols <- items.lst[["Memory"]][["MOCA"]] |> grep(pattern = "IMMT", value = TRUE)
cog_data.lst[["Memory"]] <- cog_data.lst[["Memory"]] |>
  melt(measure = cols) |>
  {
    function(DT)
    DT[
      , REGIS := sum(value), .(PTID, VISCODE2)
      ][
      , c("variable", "value") := NULL
      ]
  }() |> unique()

# Memory — MoCA: Sum of delayed recall
cols <- items.lst[["Memory"]][["MOCA"]] |> grep(pattern = "DELW", value = TRUE)
# Keep only recall w/out cue
for (col in cols) cog_data.lst[["Memory"]][
  !is.na(get(col)),
  (col) := fifelse(get(col) == 1, 1, 0)
] |> invisible()
# Sum
cog_data.lst[["Memory"]] <- cog_data.lst[["Memory"]] |>
  melt(measure = cols) |>
  {
    function(DT)
    DT[
      , DELSUM := sum(value), .(PTID, VISCODE2)
      ][
      , c("variable", "value") := NULL
      ]
  }() |> unique()

# Replace item names
items.lst[["Memory"]][["MOCA"]] <- c("REGIS", "DELSUM")

# Executive function — MoCA: Serial 7 (Total)
cols <- items.lst[["ExecFun"]][["MOCA"]] |>
  grep(pattern = "SERIAL", value = TRUE)
# Sum of Serial items
cog_data.lst[["ExecFun"]] <- cog_data.lst[["ExecFun"]] |>
  melt(measure = cols) |>
  {
    function(DT)
    DT[
      , SERIAL := sum(value), .(PTID, VISCODE2)
      ][
      , c("variable", "value") := NULL
      ]
  }() |> unique()

# Replace item names
prev_cols <- items.lst[["ExecFun"]][["MOCA"]]
items.lst[["ExecFun"]][["MOCA"]] <- c(
  prev_cols[!prev_cols %in% cols],
  "SERIAL"
)
rm(col, cols, prev_cols)


### mCFA
## Models
mods.lst    <- Map(
  function(ItemsList, FactorName) {
    items   <- paste(unlist(ItemsList), collapse = " + ")
    sprintf("%s =~ %s", FactorName, items)
  },
  items.lst,
  names(items.lst)
)

## Fits
fpath       <- here("data/rds/adni_cfa-fits_cog-domains.rds")
if (all(file.exists(fpath), !REFIT_FA)) {
  fits.lst  <- readRDS(fpath)
} else {
  fits.lst  <- Map(
    function(Model, Data) {
      # Use only most recent visit for each Subject
      fit   <- Data[order(-VISDATE)][!duplicated(PTID)] |>
      sem(
        model = Model,
        estimator = "WLSMV",
        std.lv = TRUE
      )
    },
    mods.lst,
    cog_data.lst
  )
  saveRDS(fits.lst, fpath)
}


### Apply CFA to the rest of the sessions
## Manually reviewed that the math checksout by comparing the same subjects
## using lavPredict and the computed LFactors.
## The values are not exactly the same, but the correlation is 1.0
cog_lat.dt    <- Map(
  function(Domain, Data, Fit) {
    # Unstandardized factor loadings
    loadings  <- inspect(Fit, "coef")$lambda
    # Observed variable covariance matrix
    Sigma     <- inspect(Fit, "cov.ov")
    # Factor score weights
    scr_wgts  <- solve(
      t(loadings) %*% solve(Sigma) %*% loadings
    ) %*% t(loadings) %*% solve(Sigma)
    cols      <- rownames(loadings)
    Data[
      ,
      .(
        PTID,
        VISCODE = VISCODE2,
        FACTOR = toupper(Domain),
        VAL = as.vector(as.matrix(.SD) %*% t(scr_wgts))
      ),
      .SDcols = cols
    ][!is.na(VAL)]
  },
  names(cog_data.lst),
  cog_data.lst,
  fits.lst
) |>
  rbindlist() |>
  dcast(... ~ FACTOR, value.var = "VAL")

saveRDS(cog_lat.dt, here("data/rds/adni_cfa-factors_cog-domains.rds"))

### OUTPUT
#outmds <- here('data/derivatives', paste0('mcfa-cognition_mri', 1:2, '.md'))
#for (outmd in outmds) {
  #sink(outmd)
  #parameters::parameters(mcfa_m.f) |> insight::print_md()
  #sink()
  #rm(outmd)
#}

### Save Cognition's latent factor
#outrds      <- here('data/rds/adni_cog-latent.rds')
#DT_mri[, COG_latent := lavPredict(mcfa_m.f)]
#DT_mri[, COG_latent2 := lavPredict(mcfa_m2.f)]
#cog.dt      <- DT_mri
#saveRDS(cog.dt, outrds)
