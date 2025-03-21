#!/usr/bin/env Rscript
library(here)
library(data.table)
library(Matrix)
library(Hmisc)
#library(readr)
#library(stringr)
#library(lubridate)
#library(glue)
library(gt)
library(lme4)
library(ggplot2)
library(patchwork)
library(ADNIMERGE)

ReDoPlots   <- TRUE
PrintTIFF   <- FALSE
USE_IMPUTED <- TRUE

### FUNCTIONS
here("code/functions.R") |> source()

### INPUT
fpaths      <- here(
  "data/rds",
  c(
    "adni_hc-hvr.rds",
    "adni_dxs_imputed.rds",
    "adni_age_calculated.rds",
    "adni_cog-latent.rds",
    "adni_raket_edt.rds",
    "adni_mlm-hc-hvr.rds",
    "adni_mlm-edt.rds",
    "adni_mlm-cog.rds"
  )
)

## HCvol & HVR
if (file.exists(fpaths[1])) {
  hc_hvr.dt <- readRDS(fpaths[1])
} else {
  here("code/calc_hvr_adni.R") |> source()
}

## Imputed Dxs from ADNIMERGE
if (file.exists(fpaths[2])) {
  dx.dt     <- readRDS(fpaths[2])
} else {
  here("code/impute_dx.R") |> source()
}

## AGE
## Calculated age using the date of MRI sessions
if (file.exists(fpaths[3])) {
  age.dt    <- readRDS(fpaths[3])
} else {
  here("code/calculate_age_adni.R") |> source()
}

## Cognition
# Latent variable obtained from ADASQ4, MMSE, CDRSB and validated with CFA
if (file.exists(fpaths[4])) {
  cog.dt    <- readRDS(fpaths[4])
} else {
  here("code/cfa-cog_adni_old.R") |> source()
}

### RAKET EDT
#if (file.exists(fpaths[5])) {
  #raket.dt  <- readRDS(fpaths[5])
#} else {
  #here("code/parse_raket_adni.R") |> source()
#}

## MLM
if (all(file.exists(fpaths[6:8]))) {
  fits1     <- readRDS(fpaths[6])
  fits2     <- readRDS(fpaths[7])
  fits3     <- readRDS(fpaths[8])
} else {
  here("code/model_mlm_adni.R") |> source()
}
rm(fpaths)


### Data PROCESSING
# Use imputed Dx?
if (USE_IMPUTED) dx.dt[, DX := stringr::str_remove(DX, "\\?")]
sel_dx.dt   <- dx.dt[
  ,
  .(DX = factor(DX, levels = c("CN", "MCI", "Dementia"))),
  .(PTID, EXAMDATE)]
rm(USE_IMPUTED, dx.dt)

# PTGENDER & APOE4
data(adnimerge)
setDT(adnimerge)
tiv.dt      <- adnimerge[!duplicated(PTID), .(SEX = PTGENDER, APOE4), PTID]

## HC-HVR
DT          <- tiv.dt[
  age.dt, on = "PTID"
][
  sel_dx.dt, on = .(PTID, EXAMDATE)
][
  hc_hvr.dt, on = .(PTID, EXAMDATE)
][
  ,
  .(HCv = mean(HCvol_adj), HVR = mean(HVR)),
  .(PTID, EXAMDATE, DX, AGE, SEX, APOE4)
] |> na.omit() |> setkey(PTID)

DT          <- DTclean(
  DT,
  scalevars = c("HCv", "HVR"),
  ordervars = c("APOE4", "DX"),
  centervars = "AGE",
  reference_controls = TRUE
)

# Get only AGE.c at baseline
DT          <- DT[
  .(1), on = "VIS", .(PTID, AGE.bl.c = AGE.c)
][
  DT, on = "PTID"
] |> setcolorder("AGE.bl.c", after = "AGE.bl")

# Melting
DT_l        <- DT |>
melt(
  measure = patterns("^H(Cv|VR).scl$"),
  variable.name = "HC", value.name = "VAL"
)
hcvars      <- DT_l[, levels(HC)]

# Rename vars/covars for easier formatting:
vars        <- c(
  "VAL", "TIME", "SEX", "DX", "APOE4", "AGE.bl.c", "VIS", "PTID"
)
vars_short  <- c("Y", "T", "S", "D", "A4", "A", "I", "ID")
setnames(DT_l, vars, vars_short)

### Table 1
fname <- "table-1.tex"
fpath <- here("tables")
if (!file.exists(fpath)) dir.create(fpath)
DT[, DX := factor(DX, labels = c("CH", "MCI", "AD"))]
merge(
  DT[, .SD[which.max(VIS)], by = PTID, .SDcols = "VIS"],
  DT[, .SD[which.min(VIS)], by = PTID, .SDcols = c("DX", "HCv", "HVR", "AGE")]
) |>
  melt(id = c("PTID", "DX"), variable = "VAR") |>
  suppressWarnings() |>
  na.omit() |>
  (
    \(DT) DT[
      ,
      fifelse(
        VAR == "VIS",
        sprintf("%.0f (%.2f)", mean(value), sd(value)),
        sprintf("%.2f (%.2f)", mean(value), sd(value))
      ),
      .(DX, VAR)
    ]
  )() |>
  rbind(
    DT[
      !duplicated(PTID), .N, keyby = DX
    ][
      DT[!duplicated(PTID), .N, keyby = .(DX, SEX)]
    ][
      "Female", on = "SEX",
      .(
        sprintf("%i (%.0f%%)", i.N, 100 * i.N /N),
        VAR = "SEXF"
      ),
      DX
    ],
    DT[
      !duplicated(PTID), .N, keyby = DX
    ][
      DT[!duplicated(PTID), .N, keyby = .(DX, APOE4)]
    ][
      "1", on = "APOE4",
      .(
        sprintf("%i (%.0f%%)", i.N, 100 * i.N /N),
        VAR = "APOE4_1"
      ),
      DX
    ],
    DT[
      !duplicated(PTID), .N, keyby = DX
    ][
      DT[!duplicated(PTID), .N, keyby = .(DX, APOE4)]
    ][
      "2", on = "APOE4",
      .(
        sprintf("%i (%.0f%%)", i.N, 100 * i.N /N),
        VAR = "APOE4_2"
      ),
      DX
    ],
    use.names = TRUE
  ) |>
  dcast(VAR ~ DX, value.var = "V1") |>
  (\(DT) DT[
    ,
    VAR := factor(
      VAR,
      levels = c(
        "SEXF",
        "APOE4_1",
        "APOE4_2",
        "VIS",
        "AGE",
        "HCv",
        "HVR"
      ),
      labels = c(
        "Sex (F)",
        "APOE4 (1)",
        "APOE4 (2)",
        "Visits (N)",
        "Age (bl, y)",
        "HCvol (bl, CC)",
        "HVR (bl, 0-1)"
      )
    )][order(VAR)]
  )() |>
  gt(rowname_col = "VAR", process_md = TRUE) |>
  tab_spanner(label = "Clinical Label", columns = c("CH", "MCI", "AD")) |>
  tab_options(
    latex.tbl.pos = "h",
    footnotes.multiline = FALSE
  ) |>
  cols_align("center", columns = c("CH", "MCI", "AD")) |>
  cols_label(
    CH = sprintf(
      "**CH**, N: %i", DT[!duplicated(PTID)]["CH", on = "DX", .N]
    ) |> md(),
    MCI = sprintf(
      "**MCI**, N: %i", DT[!duplicated(PTID)]["MCI", on = "DX", .N]
    ) |> md(),
    AD = sprintf(
      "**AD**, N: %i", DT[!duplicated(PTID)]["AD", on = "DX", .N]
    ) |> md()
  ) |>
  tab_footnote(
    footnote = "N (%).",
    locations = cells_stub(rows = 1:3)
  ) |>
  tab_footnote(
    footnote = "Mean (SD).",
    locations = cells_stub(rows = 4:7)
  ) |>
  gtsave(filename = fname, path = fpath)

rm(fname, fpath)

### EDT (Raket's method)
#setkey(raket.dt, PTID, EXAMDATE)
#DT_r        <- DT[
  #, .(SEX, APOE4, HCv, HVR), keyby = .(PTID, EXAMDATE)
#][raket.dt]

## Time invariant covariates
#tics.cols   <- c("PTID", "SEX", "APOE4")
#tic.dt      <- DT_r[, ..tics.cols] |>
#na.omit() |>
#unique() |>
#setkey(PTID)
#DT_r[, (tics.cols[-1]) := NULL]

## Baseline HC integrity
#bl.cols     <- c("EXAMDATE", "HCv", "HVR")
#hchvr.bl.dt <- DT_r[!is.na(HCv), .SD[which.min(AGE)],
                    #"PTID", .SDcols = bl.cols ] |>
#setnames(bl.cols, paste0(bl.cols, ".bl")) |>
#setkey(PTID)

#DT_r        <- tic.dt[hchvr.bl.dt
                      #][DT_r
                      #][EXAMDATE.bl <= EXAMDATE] |>
#unique() |>
#setcolorder(c("EXAMDATE", "MONTH", "AGE"), after = "HVR.bl")
#DT_r[, EXAMDATE.bl := NULL]

## Final CLEAN
#DT_r        <- DTclean(DT_r, scalevars = c("HCv", "HVR", "HCv.bl", "HVR.bl"),
                       #centervars = "AGE")
#DT_r        <- DT_r[VIS == 1, .(PTID, MONTH.bl = MONTH)][DT_r, on = "PTID"]
#DT_r[, MONTH := MONTH - MONTH.bl]
#DT_r[, MONTH.bl := NULL]

## Get only AGE.c at baseline
#DT_r        <- DT_r[.(1), on = "VIS", .(PTID, AGE.bl.c = AGE.c)
                    #][DT_r, on = "PTID"
                    #] |> setcolorder("AGE.bl.c", after = "AGE.bl")

## Melting
#DT_lr       <- DT_r |>
#melt(measure = patterns(VAL = "(HCv|HVR).scl", VAL.bl = "(HCv|HVR).bl.scl"),
     #variable.name = "HC")
#DT_lr[, HC := factor(HC, labels = hcvars)]

## Rename vars/covars for easier formatting:
#vars        <- c("EDT", "VAL", "VAL.bl", "AGE.bl.c", "SEX", "APOE4", "VIS", "PTID")
#vars_short  <- c("Y", "X", "Xb", "A", "S", "A4", "I", "ID")
#setnames(DT_lr, vars, vars_short)

### Cognition
bl.cols     <- c("EXAMDATE", "HCv", "HVR")
cog.dt      <- age.dt[
  cog.dt[, .(PTID, VISCODE, COG_latent)],
  on = .(PTID, VISCODE)
][
  ,
  COG_latent,
  keyby = .(PTID, EXAMDATE)
]
DT_c        <- DT[
  ,
  .(SEX, AGE, APOE4, HCv, HVR),
  keyby = .(PTID, EXAMDATE)
][cog.dt] |> setkey(PTID)

# Time invariant covariates
tics.cols   <- c("PTID", "SEX", "APOE4")
tic.dt      <- DT_c[, ..tics.cols] |>
na.omit() |>
unique() |>
setkey(PTID)
DT_c[, (tics.cols[-1]) := NULL]

# Baseline HC integrity
hchvr.bl.dt <- DT_c[
  !is.na(HCv),
  .SD[which.min(AGE)],
  "PTID",
  .SDcols = bl.cols
] |> setnames(bl.cols, paste0(bl.cols, ".bl")) |> setkey(PTID)

DT_c        <- tic.dt[hchvr.bl.dt][DT_c][EXAMDATE.bl <= EXAMDATE] |> unique()
DT_c[, EXAMDATE.bl := NULL]

# Shift Cog factor to being positive
DT_c[, COG := COG_latent + abs(min(COG_latent))]

DT_c        <- DTclean(
  DT_c,
  ordervars = "APOE4",
  scalevars = c("HCv", "HCv.bl", "HVR", "HVR.bl", "COG"),
  centervars = "AGE"
)

# AGE.c at baseline
DT_c        <- DT_c[
  .(1), on = "VIS", .(PTID, AGE.bl.c = AGE.c)
][
  DT_c, on = "PTID"
] |> setcolorder("AGE.bl.c", after = "AGE.bl")

# Melting
DT_lc       <- DT_c |> melt(
  measure = patterns(VAL = "(HCv|HVR).scl", VAL.bl = ".bl.scl"),
  variable.name = "HC"
)
DT_lc[, HC := factor(HC, labels = hcvars)]

# Rename vars/covars for easier formatting:
vars        <- c(
  "COG.scl", "VAL", "VAL.bl", "TIME", "SEX", "APOE4", "AGE.bl.c", "VIS", "PTID"
)
vars_short  <- c("Y", "X", "Xb", "T", "S", "A4", "A", "I", "ID")
setnames(DT_lc, vars, vars_short)

### Table 2
fname <- "table-2.tex"
fpath <- here("tables")
if (!file.exists(fpath)) dir.create(fpath)
merge(
  DT[
    DT_c,
    on = .(PTID, EXAMDATE),
    nomatch = NULL,
    .(VIS = .N),
    by = PTID
  ],
  DT[
    DT_c,
    on = .(PTID, EXAMDATE),
    nomatch = NULL
  ][,
    .SD[which.min(VIS)],
    by = PTID,
    .SDcols = c("DX", "HCv", "HVR", "AGE", "COG")
  ]
) |>
  melt(id = c("PTID", "DX"), variable = "VAR") |>
  suppressWarnings() |>
  na.omit() |>
  (
    \(DT) DT[
      ,
      fifelse(
        VAR == "VIS",
        sprintf("%.0f (%.2f)", mean(value), sd(value)),
        sprintf("%.2f (%.2f)", mean(value), sd(value))
      ),
      .(DX, VAR)
    ]
  )() |>
  rbind(
    DT[DT_c, on = .(PTID, EXAMDATE), nomatch = NULL][
      !duplicated(PTID), .N, keyby = DX
    ][
      DT[
        DT_c, on = .(PTID, EXAMDATE), nomatch = NULL
      ][
        !duplicated(PTID), .N, keyby = .(DX, SEX)
      ]
    ][
      "Female", on = "SEX",
      .(
        sprintf("%i (%.0f%%)", i.N, 100 * i.N /N),
        VAR = "SEXF"
      ),
      DX
    ],
    DT[DT_c, on = .(PTID, EXAMDATE), nomatch = NULL][
      !duplicated(PTID), .N, keyby = DX
    ][
      DT[
        DT_c, on = .(PTID, EXAMDATE), nomatch = NULL
      ][
        !duplicated(PTID), .N, keyby = .(DX, APOE4)
      ]
    ][
      "1", on = "APOE4",
      .(
        sprintf("%i (%.0f%%)", i.N, 100 * i.N /N),
        VAR = "APOE4_1"
      ),
      DX
    ],
    DT[DT_c, on = .(PTID, EXAMDATE), nomatch = NULL][
      !duplicated(PTID), .N, keyby = DX
    ][
      DT[
        DT_c, on = .(PTID, EXAMDATE), nomatch = NULL
      ][
        !duplicated(PTID), .N, keyby = .(DX, APOE4)
      ]
    ][
      "2", on = "APOE4",
      .(
        sprintf("%i (%.0f%%)", i.N, 100 * i.N /N),
        VAR = "APOE4_2"
      ),
      DX
    ],
    use.names = TRUE
  ) |>
  dcast(VAR ~ DX, value.var = "V1") |>
  (\(DT) DT[
    ,
    VAR := factor(
      VAR,
      levels = c(
        "SEXF",
        "APOE4_1",
        "APOE4_2",
        "VIS",
        "AGE",
        "COG",
        "HCv",
        "HVR"
      ),
      labels = c(
        "Sex (F)",
        "APOE4 (1)",
        "APOE4 (2)",
        "Visits (N)",
        "Age (bl, y)",
        "Cog. decline (lat)",
        "HCvol (bl, CC)",
        "HVR (bl, 0-1)"
      )
    )][order(VAR)]
  )() |>
  gt(rowname_col = "VAR", process_md = TRUE) |>
  tab_spanner(label = "Clinical Label", columns = c("CH", "MCI", "AD")) |>
  tab_options(
    latex.tbl.pos = "h",
    footnotes.multiline = FALSE
  ) |>
  cols_align("center", columns = c("CH", "MCI", "AD")) |>
  cols_label(
    CH = sprintf(
      "**CH**, N: %i",
      DT[
        "CH", on = "DX"
        ][
        DT_c[!duplicated(PTID)], on = .(PTID, EXAMDATE), nomatch = NULL, .N
        ]
    ) |> md(),
    MCI = sprintf(
      "**MCI**, N: %i",
      DT[
        "MCI", on = "DX"
        ][
        DT_c[!duplicated(PTID)], on = .(PTID, EXAMDATE), nomatch = NULL, .N
        ]
    ) |> md(),
    AD = sprintf(
      "**AD**, N: %i",
      DT[
        "AD", on = "DX"
        ][
        DT_c[!duplicated(PTID)], on = .(PTID, EXAMDATE), nomatch = NULL, .N
        ]
    ) |> md(),
  ) |>
  tab_footnote(
    footnote = "N (%).",
    locations = cells_stub(rows = 1:3)
  ) |>
  tab_footnote(
    footnote = "Mean (SD).",
    locations = cells_stub(rows = 4:8)
  ) |>
  gtsave(filename = fname, path = fpath)
rm(fname, fpath)


### MLM
## Chosen models from `code/model_mlm_adni.R`
## HC/HVR
hcv.mod     <- fits1[[1]][[10]]
hvr.mod     <- fits1[[2]][[10]]

DT_l["HCv.scl", on = "HC", Yp := predict(hcv.mod, DT_l["HCv.scl", on = "HC"])]
DT_l["HVR.scl", on = "HC", Yp := predict(hvr.mod, DT_l["HVR.scl", on = "HC"])]

### EDT
## Decide if use significant better model or same one for both
## Baseline HC
#edt_hcv1.mod <- fits2[[1]][[4]]
##edt_hvr1.mod <- fits2[[2]][[1]]
#edt_hvr1.mod <- fits2[[2]][[4]]


#DT_lr["HCv.scl", on = "HC",
      #Ybp := predict(edt_hcv1.mod, DT_lr["HCv.scl", on = "HC"])]
#DT_lr["HVR.scl", on = "HC",
      #Ybp := predict(edt_hvr1.mod, DT_lr["HVR.scl", on = "HC"])]

## Longitudinal HC
#edt_hcv2.mod <- fits2[[1]][[11]]
##edt_hvr2.mod <- fits2[[2]][[9]]
#edt_hvr2.mod <- fits2[[2]][[11]]

#DT_lr["HCv.scl", on = "HC",
      #Yp := predict(edt_hcv2.mod, DT_lr["HCv.scl", on = "HC"])]
#DT_lr["HVR.scl", on = "HC",
      #Yp := predict(edt_hvr2.mod, DT_lr["HVR.scl", on = "HC"])]

## Cognition
# Decide if use significant better model or same one for both
# Baseline HC
cog_hcv1.mod <- fits3[[1]][[12]]
cog_hvr1.mod <- fits3[[2]][[12]]


DT_lc[
  "HCv.scl",
  on = "HC",
  Ybp := predict(cog_hcv1.mod, DT_lc["HCv.scl", on = "HC"])
]
DT_lc[
  "HVR.scl",
  on = "HC",
  Ybp := predict(cog_hvr1.mod, DT_lc["HVR.scl", on = "HC"])
]

# Longitudinal HC
#cog_hcv2.mod <- fits3[[1]][[22]]
cog_hcv2.mod <- fits3[[1]][[24]]
cog_hvr2.mod <- fits3[[2]][[24]]

DT_lc[
  "HCv.scl",
  on = "HC",
  Yp := predict(cog_hcv2.mod, DT_lc["HCv.scl", on = "HC"])
]

DT_lc[
  "HVR.scl",
  on = "HC",
  Yp := predict(cog_hvr2.mod, DT_lc["HVR.scl", on = "HC"])
]

# Un-Z
cog.hcv.dt  <- DTunscale(
  DT_lc[!HC %like% "HVR"],
  scalevars = "Yp",
  origs = "COG",
  replace = F,
  visit.col = "I"
)

cog.hvr.dt  <- DTunscale(
  DT_lc[HC %like% "HVR"],
  scalevars = "Yp",
  origs = "COG",
  replace = F,
  visit.col = "I"
)

DT_lc["HCv.scl", on = "HC", Ypu := cog.hcv.dt$Yp.uscl]
DT_lc["HVR.scl", on = "HC", Ypu := cog.hvr.dt$Yp.uscl]

### PLOTS
## Palette
cbPalette   <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

#cbP2        <- c("#D0E8F2", "#4086AD", "#19506F", ## Gradient of blue
                 #"#F0BCC0", "#D0585E", "#7C1E23") ## Gradient of red

cbP2        <- c("#C0D8E7", "#8B8BC3", "#69186A") # APOE4
cbP3        <- c("#CCC591", "#798E87")

### HC vs HVR
f_plot1  <- sprintf(
  "plots/%s.%s", "adni_trajectory_hcv-hvr_mlm", c("png", "tiff")
) |> here()

if(
  !file.exists(f_plot1[[1]])
  | (!file.exists(f_plot1[[2]]) & PrintTIFF)
  | ReDoPlots
) {

### Factor renaming
  DTp1 <- copy(DT_l)
  DTp1[, HC := factor(HC, levels = rev(hcvars), labels = c("HVR", "Volume"))]
  DTp1[, S := factor(S, labels = c("F", "M"))]
  DTp1[, D := factor(D, labels = c("CH", "MCI", "AD"))]

  p_main <- ggplot(DTp1, aes(x = AGE, y = Y, colour = D)) +
    theme_classic(base_size = 11) +
    theme(
      text = element_text(size = 11),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    geom_line(aes(group = ID), linewidth = .25, alpha = .3) +
    #geom_smooth(aes(y = Yp), method = "lm", formula = y ~ poly(x, 2),
                #linewidth = .5, alpha = .2) +
    #geom_line(aes(group = DX, y = PRED), linewidth = .5, alpha = .2) +
    scale_colour_manual(values = cbPalette[c(4,2,7)]) +
    #scale_fill_manual(values = cbPalette[c(4,2,7)]) +
    facet_grid(rows = vars(HC)) +
    labs(subtitle = "HC vol & HVR by Age") +
    xlab("Age (years)") +
    ylab("Z-score")

  pdx <- ggplot(DTp1, aes(x = AGE, y = Y, colour = D, lty = HC)) +
    theme_classic(base_size = 11) +
    theme(
      text = element_text(size = 11),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10)
    ) +
    geom_point(size = .2, alpha = .01) +
    geom_line(
      aes(y = Yp),
      stat = "smooth",
      method = "lm",
      formula = y ~ poly(x, 2),
      linewidth = .7,
      alpha = .7
    ) +
    scale_colour_manual(values = cbPalette[c(4,2,7)]) +
    #scale_fill_manual(values = cbPalette[c(4,2,7)]) +
    facet_wrap(vars(D), ncol = 1, scales = "free") +
    labs(subtitle = "Diagnosis", colour = "Dx", lty = "HC") +
    #xlab("Age (years)") +
    ylab("Z-score")

  psex <- ggplot(DTp1, aes(x = AGE, y = Yp, colour = S, lty = HC)) +
    theme_classic(base_size = 11) +
    theme(
      text = element_text(size = 11),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10)
    ) +
    geom_line(
      stat = "smooth",
      method = "lm",
      formula = y ~ poly(x, 2),
      linewidth = .7,
      alpha = .5
    ) +
    scale_colour_manual(values = cbPalette[c(8,3)]) +
    #scale_fill_manual(values = cbPalette[c(8,3)]) +
    facet_wrap(vars(D), ncol = 1, scales = "free") +
    labs(subtitle = "Sex", colour = "Sex", lty = "HC") +
    xlab("Age (years)") #+
    #ylab("Z-score")

  pa4 <- ggplot(DTp1, aes(x = AGE, y = Yp, colour = A4, lty = HC)) +
    theme_classic(base_size = 11) +
    theme(
      text = element_text(size = 11),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10)
    ) +
    geom_line(
      stat = "smooth",
      method = "lm",
      formula = y ~ poly(x, 2),
      linewidth = .7,
      alpha = .5
    ) +
    scale_colour_manual(values = cbP2) +
    #scale_fill_manual(values = cbP2) +
    facet_wrap(vars(D), ncol = 1, scales = "free") +
    labs(subtitle = "APOE4 alleles", colour = "APOE4", lty = "HC")
    #xlab("Age (years)") +
    #ylab("Z-score")

  p1 <- p_main + (pdx + psex + pa4) +
  plot_annotation(
    title = "Neurodegeneration Trajectories",
    subtitle = paste(
      "Mixed-effect regressions:",
      "Hippocampal volume & HVR;",
      "random intercepts and slopes"
    ),
    caption = paste(
      "The cognitively healthy group was used",
      "as reference for Z scoring"
    )
  ) +
  #plot_layout(heights = c(1,1.7), guides = "collect", axes = "collect")
  plot_layout(guides = "collect", axes = "collect")

  if(!file.exists(f_plot1[[1]]) | ReDoPlots){
    ggsave(f_plot1[[1]], width = 10, height = 6, units = "in", dpi = 400)
  }

  if(PrintTIFF & (!file.exists(f_plot1[[2]]) | ReDoPlots)){
    ggsave(f_plot1[[2]], width = 7, height = 8, units = "in",
           device = "tiff", dpi = 400)
  }
}

### Cognition
f_plot3  <- sprintf(
  "plots/%s.%s", "adni_trajectory_cog_mlm", c("png", "tiff")
) |> here()

if(
  !file.exists(f_plot3[[1]])
  | (!file.exists(f_plot3[[2]]) & PrintTIFF)
  | ReDoPlots
) {
  ### Factor renaming
  DTp3 <- DT[, .(PTID, AGE, DX)][DT_lc, on = .(PTID = ID, AGE)] |>
  setnames(c("PTID", "DX"), c("ID", "D")) #|>
  #melt(measure = patterns("^X"), variable.name = "BL", value.name = "X")
  #DTp3["Xb", on = "BL", Yp := Ybp]
  #DTp3["Xb", on = "BL", Xb := X]
  #DTp3["X", on = "BL", Xl := X]
  #DTp3[, Ybp := NULL]
  #DTp3[, BL := factor(BL, levels = c("Xb", "X"),
                      #labels = c("Baseline only", "All visits"))]
  DTp3[, D := factor(D, labels = c("CH", "MCI", "AD"))]
  DTp3[, HC := factor(HC, levels = rev(hcvars), labels = c("HVR", "Volume"))]
  DTp3[, S := factor(S, labels = c("F", "M"))]


  p_main <- DTp3["HVR", on = "HC"] |>
  ggplot(aes(x = AGE, y = COG, colour = D)) +
    theme_classic(base_size = 11) +
    theme(
      text = element_text(size = 11),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      legend.position = "none"
    ) +
    geom_line(aes(group = ID), linewidth = .35, alpha = .25) +
    geom_hline(yintercept = 0, lty = "dashed", colour = cbPalette[1]) +
    scale_colour_manual(values = cbPalette[c(4,2,7)]) +
    labs(subtitle = "Individual trajectories") +
    #xlab("Age (years)") +
    ylab("Decline (observed)")

  pgroup <- DTp3[AGE > 57 & AGE < 93]["HVR", on = "HC"] |>
  ggplot(aes(x = AGE, y = COG)) +
    theme_classic(base_size = 11) +
    theme(
      text = element_text(size = 11),
      #axis.title.y = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      legend.position = "right"
    ) +
    scale_fill_manual(values = cbPalette[c(4,2,7)]) +
    scale_colour_manual(values = cbPalette[c(4,2,7)]) +
    geom_hline(yintercept = 0, lty = "dashed", colour = cbPalette[1]) +
    geom_smooth(
      aes(colour = D, fill = D),
      method = "lm",
      formula = y ~ poly(x, 2),
      linewidth = .5,
      se = T,
      alpha = .2
    ) +
    geom_smooth(
      colour = cbPalette[1],
      fill = cbPalette[1],
      method = "lm",
      formula = y ~ poly(x, 2),
      linewidth = .7,
      se = T,
      alpha = .2
    ) +
    labs(subtitle = "Group trajectories", colour = "Dx", fill = "Dx") +
    ylab("Decline (observed)") +
    xlab("Age (years)")


  #phc <- ggplot(aes(x = COG, y = X, colour = HC, fill = HC, lty = HC)) +
  phc <-  DTp3["HVR", on = "HC"] |>
    ggplot(aes(x = COG, y = X, colour = HC, fill = HC)) +
    theme_classic(base_size = 11) +
    theme(
      text = element_text(size = 11),
      #axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      legend.position = "none"
    ) +
    geom_point(size = .5, shape = 21, fill = "transparent", alpha = .15) +
    geom_smooth(
      aes(x = Ypu),
      linewidth = .7,
      alpha = .1,
      method = "lm",
      formula = y ~ poly(x, 2)
    ) +
    scale_colour_manual(values = cbPalette[1]) +
    scale_fill_manual(values = cbPalette[1]) +
    labs(subtitle = "Cognition & HC",
         x = "Decline (predicted)", y = "HC: Z-score") +
    xlim(0, 5)

  psex <- DTp3["HVR", on = "HC"] |>
    ggplot(aes(x = Ypu, y = X, colour = S, lty = HC)) +
    theme_classic(base_size = 11) +
    theme(
      text = element_text(size = 11),
      #axis.title.y = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10)
    ) +
    geom_line(
      stat = "smooth",
      method = "lm",
      formula = y ~ poly(x, 2),
      linewidth = .7,
      alpha = .5
    ) +
    scale_colour_manual(values = cbPalette[c(8,3)]) +
    #scale_fill_manual(values = cbPalette[c(8,3)]) +
    labs(subtitle = "Cognition & Sex", colour = "Sex",
         x = "Decline (predicted)", y = "HC: Z-score") +
    guides(lty = "none") +
    xlim(0, 5)
    #xlab("Decline (predicted)") #+
    #ylab("Z-score")


  pa4 <- DTp3["HVR", on = "HC"] |>
    ggplot(aes(x = Ypu, y = X, colour = A4, lty = HC)) +
    theme_classic(base_size = 11) +
    theme(
      text = element_text(size = 11),
      #axis.title.x = element_blank(),
      #axis.title.y = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10)
    ) +
    geom_line(
      stat = "smooth",
      method = "lm",
      formula = y ~ poly(x, 2),
      linewidth = .7,
      alpha = .5
    ) +
    scale_colour_manual(values = cbP2) +
    #scale_fill_manual(values = cbP2) +
    #facet_wrap(vars(D), ncol = 1) +
    guides(lty = "none") +
    labs(subtitle = "Cognition & APOE4", colour = "APOE4",
         x = "Decline (predicted)", y = "HC: Z-score") +
    xlim(0, 5)

  p3 <- ((p_main / pgroup) | (phc / psex / pa4)) +
    plot_annotation(
      title = "Cognitive Decline Trajectories",
      subtitle = paste(
        "Mixed-effect regressions:",
        "HVR;",
        "random intercepts and slopes"
      ),
      caption = paste(
        "Cognitive decline was measured",
        "as a latent variable obtained from:",
        "ADASQ4, MMSE & CDRSB"
      )
    ) +
    #plot_layout(heights = c(1, .6, .7), guides = "collect", axes = "collect")
    plot_layout(widths = c(1, .5), guides = "collect", axes = "collect")

  if(!file.exists(f_plot3[[1]]) | ReDoPlots){
    ggsave(f_plot3[[1]], width = 10, height = 6, units = "in", dpi = 400)
  }

  if(PrintTIFF & (!file.exists(f_plot3[[2]]) | ReDoPlots)){
    ggsave(f_plot3[[2]], width = 7, height = 8, units = "in",
           device = "tiff", dpi = 400)
  }
}
