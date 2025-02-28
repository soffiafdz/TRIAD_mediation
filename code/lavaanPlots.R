## Path plots
fname     <- here("plots/mediation_paths_hvr_ug_full.pdf")
if (!file.exists(fname) & print_plots) {
  lavaanPlot2(model = mod_hvr, labels = c(SEX_n = "SEX"),
              graph_options = list(rankdir = "LR"),
              node_options = list(shape = "box"),
              edge_options = list(color = "grey"),
              coef_labels = T, stand = T, stars = "regress") |>
       embed_plot_pdf(fname)
}

fname     <- here("plots/mediation_paths_hvr_ug_eskeleton.pdf")
if (!file.exists(fname) & print_plots) {
  lavaanPlot2(model = mod_hvr, labels = c(SEX_n = "SEX"),
              graph_options = list(rankdir = "LR"),
              node_options = list(shape = "box"),
              edge_options = list(color = "grey"),
              coef_labels = F) |>
       embed_plot_pdf(fname)
}

fname     <- here("plots/mediation_paths_hvr_ug_sign.pdf")
if (!file.exists(fname) & print_plots) {
  coefs   <- extract_coefs(mod_hvr, stand = TRUE) |>
    as.data.table() |>
    {\(x) x[p_val < 0.05]}()
  if (coefs[, .N] == 0) next
  ndf     <- create_nodes(coefs, c(SEX_n = "SEX"), list(shape = "box"))
  edf     <- create_edges(coefs, ndf, list(color = "grey"),
                          coef_labels = TRUE, stars = "regress")
  dot     <- convert_graph(ndf, edf, list(rankdir = "LR"))
  lavaanPlot2(gr_viz = dot) |>
    embed_plot_pdf(fname)
  #rm(coefs, ndf, edf, dot)
}

## HCv full models
fname       <- here("data/rds/med-mod_hcv_raket_ug.rds")
if (!file.exists(fname) | refit_mods) {
  mod_hcv <-
    sem(hcv.mod,
        data = DT,
        estimator = "ML",
        se = "bootstrap",
        bootstrap = 1000)

  write_rds(mod_hcv, fname)
} else {
  mod_hcv   <- read_rds(fname)
}
rm(fname)

## Path plots
fname     <- here("plots/mediation_paths_hcv_ug_full.pdf")
if (!file.exists(fname) & print_plots) {
  lavaanPlot2(model = mod_hcv, labels = c(SEX_n = "SEX"),
              graph_options = list(rankdir = "LR"),
              node_options = list(shape = "box"),
              edge_options = list(color = "grey"),
              coef_labels = T, stand = T, stars = "regress") |>
       embed_plot_pdf(fname)
}

fname     <- here("plots/mediation_paths_hcv_ug_sign.pdf")
if (!file.exists(fname) & print_plots) {
  coefs   <- extract_coefs(mod_hcv, stand = TRUE) |>
    as.data.table() |>
    {\(x) x[p_val < 0.05]}()
  if (coefs[, .N] == 0) next
  ndf     <- create_nodes(coefs, c(SEX_n = "SEX"), list(shape = "box"))
  edf     <- create_edges(coefs, ndf, list(color = "grey"),
                          coef_labels = TRUE, stars = "regress")
  dot     <- convert_graph(ndf, edf, list(rankdir = "LR"))
  lavaanPlot2(gr_viz = dot) |>
    embed_plot_pdf(fname)
  #rm(coefs, ndf, edf, dot)
}

