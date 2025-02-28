#!/usr/bin/env Rscript

library(here)
library(data.table)

### INPUT
fpaths <- c("NEUROBAT", "ADAS_ADNIGO23", "MMSE", "MOCA") |>
    sprintf(fmt = "data/adni/cognitive_tests/%s_Jan2025.csv") |>
    here()

cogdata.lst <- fpaths |>
  lapply(\(fpath) {
    if (!file.exists(fpath)) {
      sprintf(
        "File: %s is required but could not be found",
        basename(fpath)
      ) |> stop()
    }

    fpath |>
      fread(drop = "VISCODE") |>
      setnames("VISCODE2", "VISCODE") |>
      setkey(PTID, VISCODE)
  }) |>
  setattr(
    "names",
    fpaths |>
      lapply(basename) |>
      lapply(sub, pattern = "_.*", replacement = "")
  )

rm(fpaths)

## Load ADNIMERGE
library(ADNIMERGE)
data(adnimerge)
setDT(adnimerge)
DT <- adnimerge[, .(PTID, VISCODE)] |> setkey(PTID, VISCODE)
rm(adnimerge)


### Cognitive DOMAINS
## Specify the items necessary to construct 4 latent factors
## for Memory, Executive function, Language & Visuospatial function
## from NEUROBAT, ADAS13, MMSE & MoCA
## on the ADNI data

## List for cognitive domains
cogdomains.lst   <- list(
  Memory = list(
    NEUROBAT = c(
      # WMS-R: Logical Memory - immediate & recall (Continuous)
      "LIMMTOTAL", "LDELTOTAL",
      # Rey-AVLT:  Trials 1:6 (Continuous)
      "AVTOT1", "AVTOT2", "AVTOT3", "AVTOT4", "AVTOT5", "AVTOT6",
      # Rey-AVLT:  Trial B, 30 min delay & recognition (Continuous)
      "AVTOTB", "AVDEL30MIN", "AVDELTOT"
    ),
    ADAS = c(
      "Q1SCORE", # Word recall (Continuous)
      "Q4SCORE", # Delayed word recall (Continuous - 11 levels)
      "Q7SCORE", # Orientation (Continuous - 9 levels)
      "Q8SCORE"  # Word recognition (Continuous - 13 levels)
    ),
    MMSE = c( ## All are binary
      # Orientation :: Sum all of these
      "MMDATE", "MMYEAR", "MMMONTH", "MMDAY", "MMSEASON",
      "MMHOSPIT", "MMFLOOR", "MMCITY", "MMAREA", "MMSTATE",
      # Word (Ball): immediate/delayed recall
      # These are binary, but only two, so sum is worthless
      "WORD1", "WORD1DL"
    ),
    MOCA = c(
      ## Sums and Averages are continuous
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
      # Clock (Binary) :: Sum all of these
      "CLOCKCIRC", "CLOCKSYM", "CLOCKNUM", "CLOCKHAND", "CLOCKTIME",
      # WAIS-R: Digit span backward
      #"DSPANBAC",
      ## WMS-R: Digit symbol & Digit span forward
      #"DIGITSCOR", "DSPANFOR"
      # Trails A/B Time (Continuous)
      "TRAASCOR", "TRABSCOR"
    ),
    ADAS = "Q13SCORE",    # Number cancellation task (Continuous? - 6 levels)
    #MMSE = "WORLDSCORE",  # Spell `world` backwards
    MOCA = c(
      # Abstraction (Binary)
      "ABSMEAS", "ABSTRAN",
      # Trails (Binary)
      "TRAILS",
      # Digits backward/forward (Binary)
      "DIGBACK", "DIGFOR",
      ## This needs to be cleaned ##
      ## Sum, but 5 levels. Discrete??
      # Serial 7 (Total)
      "SERIAL1", "SERIAL2", "SERIAL3", "SERIAL4", "SERIAL5",
      ## This needs to be cleaned ##
      # Letters/tapping (Errors) (Continuous)
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
      "Q2SCORE",  # Commands (Discrete? - 5 levels)
      "Q5SCORE",  # Naming (Continuous? - 6 levels)
      "Q6SCORE"   # Ideational-Praxis (Continuous? - 6 levels)
    ),
    MMSE = c(     ## All are binary
    "MMREPEAT",   # Repeat after instructor
    "MMHAND",     # Take paper
    "MMFOLD",     # Fold paper
    "MMONFLR",    # Place paper on floor
    "MMREAD",     # Read paper
    "MMWRITE"     # Write a sentence
    ),
    MOCA = c(
      # Naming (Binary)
      "CAMEL", "LION", "RHINO",
      # Repeat sentence (Binary)
      "REPEAT1", "REPEAT2",
      # Fluency-F (Continuous)
      "FFLUENCY"
    )
  ),
  VisospatialFun = list(
    NEUROBAT = c(
      # Clock (Binary)
      "COPYCIRC", "COPYSYM", "COPYNUM", "COPYTIME"
    ),
    ADAS = "Q3SCORE", # Constructional Praxis (Discrete? - 5 levels)
    MMSE = "MMDRAW"   # Interlocking pentagons (Binary)
  )
)

# List of ordinal/binary data
## List for cognitive domains
ordinal_items.lst   <- list(
  MMSE = c( ## All are binary
    # Orientation
    "MMDATE", "MMYEAR", "MMMONTH", "MMDAY", "MMSEASON",
    "MMHOSPIT", "MMFLOOR", "MMCITY", "MMAREA", "MMSTATE",
    # Word (Ball): immediate/delayed recall
    "WORD1", "WORD1DL",
    "MMREPEAT",   # Repeat after instructor
    "MMHAND",     # Take paper
    "MMFOLD",     # Fold paper
    "MMONFLR",    # Place paper on floor
    "MMREAD",     # Read paper
    "MMWRITE",    # Write a sentence
    "MMDRAW"      # Interlocking pentagons
  ),
  NEUROBAT = c(
    # Clock (Binary)
    "CLOCKCIRC", "CLOCKSYM", "CLOCKNUM", "CLOCKHAND", "CLOCKTIME",
    "COPYCIRC", "COPYSYM", "COPYNUM", "COPYTIME"
  ),
  MOCA = c(
    # Abstraction (Binary)
    "ABSMEAS", "ABSTRAN",
    # Trails (Binary)
    "TRAILS",
    # Digits backward/forward (Binary)
    "DIGBACK", "DIGFOR",
    # Naming (Binary)
    "CAMEL", "LION", "RHINO",
    # Repeat sentence (Binary)
    "REPEAT1", "REPEAT2"
  )#,
  #ADAS = c(
    #"Q2SCORE",  # Commands (Discrete? - 5 levels)
    #"Q3SCORE",  # Constructional Praxis (Discrete? - 5 levels)
    #"Q5SCORE",  # Naming (Continuous? - 6 levels)
    #"Q6SCORE",  # Ideational-Praxis (Continuous? - 6 levels)
    #"Q13SCORE"  # Number cancellation task (Continuous? - 6 levels)
  #)
)

## Create a list for all cognitive domains' items
testnames <- cogdomains.lst |> lapply(names) |> unlist() |> unique()
items.lst <- lapply(testnames, \(testname) {
  items <- lapply(cogdomains.lst, \(domain) {domain[[testname]]})
  items |> unlist(use.names = FALSE) |> unique() |> sort()
}) |> setattr("names", testnames)



## Aggregate cognitive scores to single list
cog.lst <- list()
for (test in testnames) {
  items   <- c("PTID", "VISCODE", "VISDATE", items.lst[[test]])
  cog.lst[[test]] <- cogdata.lst[[test]][
    , ..items
  ][
    DT, on = .(PTID, VISCODE), nomatch = NULL
  ][
    , let(
      VISCODE = NULL,
      VISDATE = lubridate::ymd(VISDATE)
    )
  ]

  setnames(
    cog.lst[[test]],
    items[-(1:2)],
    sprintf("%s_%s", test, tolower(items[-(1:2)]))
  )
  rm(items, test)
}

rm(testnames, items.lst, cogdata.lst)

## Cognitive Data cleaning
# Memory — MoCA: Sum of trials of registration
cols <- cog.lst[["MOCA"]] |>
  names() |>
  grep(pattern = "MOCA_immt", value = TRUE)

cog.lst[["MOCA"]] <- cog.lst[["MOCA"]] |>
  melt(measure = cols) |>
  {\(DT) DT[
    , MOCA_regis := sum(value), .(PTID, MOCA_visdate)
  ][
    , c("variable", "value") := NULL
  ]}() |>
    unique()

# Memory — MoCA: Sum of delayed recall
cols <- cog.lst[["MOCA"]] |>
  names() |>
  grep(pattern = "MOCA_delw", value = TRUE)

cog.lst[["MOCA"]] <- cog.lst[["MOCA"]][
  ,
  # Keep only recall w/out cue
  (cols) := lapply(.SD, \(x) fifelse(x == 1, 1, 0)),
  .SDcols = cols
] |>
  melt(measure = cols) |>
  { \(DT) DT[
    , MOCA_delsum := sum(value), .(PTID, MOCA_visdate)
  ][
    , c("variable", "value") := NULL
  ]}() |>
    unique()

# Executive function — MoCA: Serial 7 (Total)
cols <- cog.lst[["MOCA"]] |>
  names() |>
  grep(pattern = "MOCA_ser", value = T) |>
  sort()

cog.lst[["MOCA"]] <- cog.lst[["MOCA"]] |>
  melt(measure = cols) |>
  {\(DT) DT[
    , MOCA_serial := sum(value), .(PTID, MOCA_visdate)
  ][
    , c("variable", "value") := NULL
  ]}() |>
    unique()

rm(cols)

## Reorder columns for MoCA
cog.lst[["MOCA"]] <- cog.lst[["MOCA"]] |>
  names() |>
  grep(pattern = "MOCA_", value = TRUE) |>
  sort() |>
  grep(pattern = "visdate", invert = TRUE) |>
  setcolorder(x = cog.lst[["MOCA"]], after = "MOCA_visdate")

## Apply changes to the cogdomains list
cogdomains.lst[["Memory"]][["MOCA"]] <- c("REGIS", "DELSUM")
cogdomains.lst[["ExecFun"]][["MOCA"]] <- c(
  grep(
    "SERIAL",
    cogdomains.lst[["ExecFun"]][["MOCA"]],
    value = TRUE,
    invert = TRUE
  ),
  "SERIAL"
)

## Unlist cogdomains.lst
cogdomains.lst <- lapply(
  cogdomains.lst,
  \(domain) {
    Map(
      \(name, items) sprintf("%s_%s", name, tolower(items)),
      names(domain),
      domain
    ) |>
    unlist(use.names = FALSE) |>
    sort()
  }
)

### OUTPUT
fpaths <- c("data", "domains_items", "items-ordinal") |>
  sprintf(fmt = "data/rds/adni_cognitive-%s.rds") |>
  here()

saveRDS(cog.lst, fpaths[1])
saveRDS(cogdomains.lst, fpaths[2])
saveRDS(ordinal_items.lst, fpaths[3])
rm(fpaths)
