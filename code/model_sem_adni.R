#!/usr/bin/env Rscript

library(here)
library(data.table)
library(ADNIMERGE)
library(labelled)
library(lubridate)
library(lavaan)
library(lavaanPlot)
#library(lme4)
#library(stringr)
#library(stargazer)
#library(progress)
#library(parameters)

REFITMODELS <- F
PRINTPLOTS <- T

### INPUT
fpaths <- list(
  CSV = here("data/adni_wmh-vol.csv"),
  #COGTST = c("NEUROBAT", "ADAS_ADNIGO23", "MMSE", "MOCA") |>
    #sprintf(fmt = "data/adni/cognitive_tests/%s_Jan2025.csv") |>
    #here(),
  RDS = c(
  "adni_hc-hvr",
  "adni_dxs_imputed",
  "adni_age_calculated",
  "ucb_pet-amy",
  "ucb_pet-tau",
  "adni_cognitive-data",
  "adni_cognitive-domains_items",
  "adni_cognitive-items-ordinal"
  #"../adni_wmh-vol",
  #"adni_cog-latent",
  #"adni_biomarkers-latent"
  ) |> sprintf(fmt = "data/rds/%s.rds") |> here(),
  SRC = c(
    "calc_hvr_adni",
    "impute_dx",
    "calculate_age_adni",
    "parse_pet_adni",
    "spec_cogdomains"
    #"cfa-cog_adni"
    #"cfa-cog_adni",
    #"cfa-biomarkers_adni"
  ) |> sprintf(fmt = "code/%s.R") |> here()
)

## HCvol & HVR
data.lst <- list()
if (file.exists(fpaths[["RDS"]][1])) {
  #hc_hvr.dt <- fpaths[["RDS"]][1] |> readRDS()
  data.lst[["HC_HVR"]] <- fpaths[["RDS"]][1] |> readRDS()
} else {
  fpaths[["SRC"]][1] |> source()
  data.lst[["HC_HVR"]] <- hc_hvr.dt
}

## Imputed Dxs from ADNIMERGE
if (file.exists(fpaths[["RDS"]][2])) {
  #dx.dt     <- fpaths[["RDS"]][2] |> readRDS()
  data.lst[["DX"]] <- fpaths[["RDS"]][2] |> readRDS()
} else {
  fpaths[["SRC"]][2] |> source()
  data.lst[["DX"]] <- dx.dt
}

## AGE
## Calculated age using the date of MRI sessions
if (file.exists(fpaths[["RDS"]][3])) {
  #age.dt    <- fpaths[["RDS"]][3] |> readRDS()
  data.lst[["AGE"]] <- fpaths[["RDS"]][3] |> readRDS()
} else {
  fpaths[["SRC"]][3] |> source()
  data.lst[["AGE"]] <- age.dt
}

## PET
if (all(file.exists(fpaths[["RDS"]][4:5]))) {
  data.lst[["PET"]] <- lapply(fpaths[["RDS"]][4:5], readRDS) |>
    setattr("names", c("AMY", "TAU"))
  #ucb_amy.dt <- fpaths[["RDS"]][4] |> readRDS()
  #ucb_tau.dt <- fpaths[["RDS"]][5] |> readRDS()
} else {
  fpaths[["SRC"]][4] |> source()
  data.lst[["PET"]] <- list(AMY = ucb_amy.dt, TAU = ucb_tau.dt)
}

## WMH
if (file.exists(fpaths[["CSV"]])) {
  # TODO: Check that gsub works
  data.lst[["WMH"]] <- fpaths[["CSV"]] |> fread(key = c("PTID", "EXAMDATE"))
} else {
  sprintf("File: %s is required but could not be found.", fpaths[6]) |> stop()
}

## Cognition
if (all(file.exists(fpaths[["RDS"]][6:8]))) {
  data.lst[["COG"]] <- lapply(fpaths[["RDS"]][6:8], readRDS) |>
    setattr("names", c("SCRS", "DOMS", "ORD"))
} else {
  fpaths[["SRC"]][5] |> source()
  data.lst[["COG"]] <- list(
    SCRS = cog.lst,
    DOMS = cogdomains.lst,
    ORD = ordinal_items.lst
  )
}

rm(fpaths)

### Data PROCESSING
## Get PTGENDER & APOE4 from ADNIMERGE
data(adnimerge)
setDT(adnimerge)
data.lst[["COVARS"]] <- adnimerge |>
  (
    \(DT)
    DT <- DT[
      !duplicated(PTID),
      .(
        PTID,
        SEX = fifelse(PTGENDER %like% "F", 1, 0),
        APOE4,
        EDUC = PTEDUCAT
      )
    ][, SX_A4 := SEX * APOE4]
  )() |>
  merge(data.lst[["AGE"]], by = "PTID") |>
  setkey(PTID, EXAMDATE)
rm(adnimerge)

## Base Data.table
# MERGE and average by side
data.dt <- data.lst[["COVARS"]][
  data.lst[["HC_HVR"]][
    ,
    .(HCv = mean(HCvol_adj), HVRi = 1 - mean(HVR)),
    keyby = .(PTID, EXAMDATE)
  ],
  nomatch = NULL
#][
  #data.lst[["WMH"]],   ## Am I exploring WMH???
  #nomatch = NULL
][
  data.lst[["PET"]][["AMY"]][
    abs(DIFF_amy_mri) < 365,
    .(AMY_diff = DIFF_amy_mri, AMY_tracer = TRACER, AMY = CENTILOIDS),
    keyby = .(PTID, EXAMDATE = DATE_MRI)
  ],
  nomatch = NULL
][
  data.lst[["PET"]][["TAU"]][
    abs(DIFF_tau_mri) < 365,
    .(TAU_diff = DIFF_tau_mri, TAU = META_TEMPORAL_SUVR),
    keyby = .(PTID, EXAMDATE = DATE_MRI)
  ],
  nomatch = NULL
]

## Cognitive data
data.lst[["COG"]][["DATES"]] <- lapply(
  data.lst[["COG"]][["SCRS"]],
  \(DT) {
    cols <- grep("PTID|visdate", names(DT), value = TRUE)
    DT[
      , ..cols
    ][
      data.dt, on = "PTID", allow.cartesian = TRUE, nomatch = NULL
    ][
      , .SD[which.min(abs(get(cols[2]) - EXAMDATE))], .(PTID, VISCODE)
    ][
      ,
      .(VISDATE = get(cols[2]), DIFF = get(cols[2]) - EXAMDATE),
      .(PTID, VISCODE)
    ][
      , .SD[which.min(abs(DIFF))], .(PTID, VISDATE)
    ][
    abs(DIFF) < 365
    ] |>
    setcolorder("VISCODE", after = "PTID") |>
    setnames("VISDATE", cols[2]) |>
    setnames("DIFF", sub("visdate", "diff", cols[2]))
  }
)

data.lst[["COG"]][["SCRS"]] <- Map(
  \(scores, dates){
    cols <- names(scores)[1:2]
    setkeyv(scores, cols)
    setkeyv(dates, cols)
    dates[scores, nomatch = NULL] |> setkey(PTID, VISCODE)
  },
  data.lst[["COG"]][["SCRS"]],
  data.lst[["COG"]][["DATES"]]
)

ord_items <- data.lst[["COG"]][["ORD"]] |>
  (\(ord_items.lst) {
    testnames <- names(ord_items.lst)
    Map(
      \(test, items) sprintf("%s_%s", test, tolower(items)),
      testnames,
      ord_items.lst
    )
  })() |>
  unlist(use.names = FALSE)

clean_data.lst <- lapply(
  data.lst[["COG"]][["DOMS"]],
  \(items) {
    DT <- copy(data.dt)
    setkey(DT, PTID, VISCODE)
    for(subDT in data.lst[["COG"]][["SCRS"]]) {
      diffcol <- grep("diff", names(subDT), value = TRUE)
      cols <- c(
        "PTID", "VISCODE", diffcol, names(subDT)[names(subDT) %in% items]
      )
      DT <- merge(DT, subDT[, ..cols], all.x = TRUE) |>
        unique() |>
        setnames(diffcol, sub("_", "__", diffcol))
    }
    DT[
      ,
      COG_avg_diff := rowMeans(sapply(.SD, as.numeric)),
      .SDcols = patterns("__")
    ]
    DT[, (grep("__", names(DT), value = TRUE)) := NULL]
    setcolorder(DT, "COG_avg_diff", after = "TAU")
    #na.omit(DT)
    DT <- na.omit(DT)
    items_ord <- ord_items[ord_items %in% items]
    DT[, (items_ord) := lapply(.SD, ordered), .SDcols = items_ord]
  }
)

### All in a single data.table
#setkey(data.dt, PTID, VISCODE)
#for (DT in data.lst[["COG"]][["SCRS"]]) {
  #data.dt <- merge(data.dt, DT, all.x = TRUE) |> unique()
#}
#rm(DT)

### MODEL specification
### Ordinal data in a Multilevel SEM is not supported in Lavaan.
### I will specify two implementations of the models:
### 1) will respect the ordinal/binary nature of the indicators
### 2) will implement the between/within models.

## Structural Model (preliminary)
str_mod <- '
AMY  ~ a1*AGE + a2*SEX + a3*APOE4 + a4*SX_A4 + a5*AMY_diff + a6*AMY_tracer
TAU  ~ a*AMY + t1*AGE + t2*SEX + t3*APOE4 + t4*SX_A4 + t5*TAU_diff
HVRi ~ b*AMY + c*TAU + h1*AGE + h2*SEX + h3*APOE4 + h4*SX_A4
%s   ~ d*AMY + e*TAU + f*HVRi + m1*AGE + m2*EDUC + m3*SEX + m4*COG_avg_diff
'

## Mediation derived effects
med_mod <- '
Total   := d + (a*e) + (b*f) + (a*c*f)
deAMY   := d
ieTAU   := a * e
ieHVR   := b * f
ieSerial:= a * c * f
ieTotal := (a*e) + (b*f) + (a*c*f)
propAMY := deAMY / Total
propTAU := ieTAU / Total
propHVR := ieHVR / Total
propSerial := ieSerial / Total
propTotal  := ieTotal / Total
'

## Measurement model
models.lst <- Map(
  \(domain, items){
    list(
      MSR = sprintf("%s =~ %s", domain, paste(items, collapse = " + ")),
      STR = sprintf(str_mod, domain),
      MED = med_mod
    )
  },
  names(data.lst[["COG"]][["DOMS"]]),
  data.lst[["COG"]][["DOMS"]]
) |>
  lapply(paste, collapse = "\n")

#lapply(data.lst[["COG"]][["DOMS"]],
  #\(sublist)
  #paste(
    #domain,
    #"=~",
    #paste(
      #data.lst[["COG"]][["DOMS"]][[domain]],
      #collapse = " + "
    #)
  #)
#)

## Structural model
#structural_mod <- paste(
  #"Cognition",
  #"=~",
  #paste(cog_domains, collapse = " + ")
#)

### Model fitting
fits.lst <- Map(
  \(model, DT){
    sem(
      model = model,
      data = DT,
      estimator = "WLSMV",
      #cluster = "PTID",
      #control = list(iter.max = 20000)
      std.lv = T
    )
  },
  models.lst[1],
  clean_data.lst[1]
)

