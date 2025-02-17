#!/usr/bin/env Rscript

library(data.table)

### FUNCTIONS
## Data CLEANING
DTclean     <- function(
  DT,
  scalevars = NULL,
  ordervars = NULL,
  centervars = NULL,
  time_from_age = T,
  scale_baseline = T,
  reference_controls = FALSE,
  ref_column = "DX",
  ref_label = "CN",
  center_baseline = T
) {
  if (time_from_age) {
    DTw <- copy(DT)
    DTw[, c("AGE.bl", "TIME") := NULL] |> suppressWarnings()
    DTw[order(AGE), VIS := rowid(PTID)]
    setcolorder(DTw, "VIS", after = "PTID")
    DTw <- DTw[DTw[.(1), on = "VIS", .(AGE.bl = AGE), "PTID"], on = "PTID"]
    DTw[, TIME := AGE - AGE.bl]
    setcolorder(DTw, c("AGE.bl", "TIME"), before = "AGE")
  } else {
    scale_baseline <- center_baseline <- FALSE
  }
  if (!is.null(ordervars)) {
    DTw[, (ordervars) := lapply(.SD, as.ordered), .SDcols = ordervars]
  }
  if (!is.null(centervars)) {
    if (center_baseline) {
      DTp <- DTw[
        .(1),
        on = "VIS",
        lapply(.SD, mean, na.rm = T),
        .SDcols = centervars
      ]
      DTw[
        ,
        (paste0(centervars, ".c")) := lapply(
          centervars, \(x) {(get(x) - DTp[[x]][1])}
        )
      ]
    } else {
      DTw[
        ,
        (paste0(centervars, ".c")) := lapply( .SD, scale, scale = F),
        .SDcols = centervars
      ]
    }
  }
  if (!is.null(scalevars)) {
    ##TODO: Find out why this doesn't work.
    ## Have to use a vector of characters.
    ##if (is.integer(scalevars)) scalevars <- names(DTw)[scalevars]
    #browser()
    DTr <- if (scale_baseline) DTw[.(1), on = "VIS"] else DTw
    DTr <- if (reference_controls) DTr[get(ref_column) == ref_label] else DTr
    DTp <- rbind(
      DTr[, lapply(.SD, mean, na.rm = TRUE), .SDcols = scalevars],
      DTr[, lapply(.SD, sd, na.rm = TRUE), .SDcols = scalevars]
    )
    DTw[
      ,
      (paste0(scalevars, ".scl")) := lapply(
        scalevars, \(x) {(get(x) - DTp[[x]][1]) / DTp[[x]][2]}
      )
    ]
  }
  return(DTw)
}

## Bring back from Z scores
DTunscale <- function(
  DT,
  scalevars,
  origs,
  scale_baseline = TRUE,
  replace = TRUE,
  newnames = NULL,
  visit.col = "ID"
) {
  DTw <- copy(DT)
  if (scale_baseline) {
    DTp <- rbind(
      DTw[.(1), on = visit.col, lapply(.SD, mean), .SDcols = origs],
      DTw[.(1), on = visit.col, lapply(.SD, sd), .SDcols = origs]
    )
  } else {
    DTp <- rbind(
      DTw[, lapply(.SD, mean), .SDcols = origs],
      DTw[, lapply(.SD, sd), .SDcols = origs]
    )
  }
  setnames(DTp, scalevars)
  if (replace) {
    DTw[
      ,
      (scalevars) := lapply(
        scalevars, \(x) {(get(x) * DTp[[x]][2]) + DTp[[x]][1]}
      )
    ]
  } else {
    newnames <- if (!is.null(newnames)) {
      newnames
    } else if (any(stringr::str_detect(scalevars, "scl"))) {
      sub("scl", "uscl", scalevars)
    } else {
      paste0(scalevars, ".uscl")
    }
    DTw[
      ,
      (newnames) := lapply(
        scalevars, \(x) {(get(x) * DTp[[x]][2]) + DTp[[x]][1]}
      )
    ]
  }
}

## build regression FORMULAS
build_formulas <- function(
  Y,
  X,
  interactionvars = NULL,
  notinteractionvars = NULL,
  random_effects = c("none", "intercepts", "slopes"),
  id = NULL,
  slopevars = NULL,
  itervars = NULL,
  quadraticvars = NULL,
  quadratic_interactions = FALSE,
  skip_items = NULL
) {
  # Quadratic terms
  if (is.null(quadraticvars)) {
    quadratic_interactions <- FALSE
  } else {
    quad_terms <- sprintf("I(%s^2)", quadraticvars)
    if (quadratic_interactions) {
      quad_inters <- lapply(
        quad_terms, paste,
        X[!X %in% quadraticvars],
        sep = ":"
      ) |> unlist()
    }
  }

  # Iteration of terms
  if (!is.null(itervars)) {
    iter_terms <- itervars
    if (length(itervars) > 1) {
      for (i in 2:length(itervars)) {
        iters <- combn(itervars, i) |>
        t() |> apply(1, \(x) {paste(x, collapse = "+")})
        iter_terms <- c(iter_terms, iters)
      }
    }
    if (!is.null(skip_items)) skip_items <- skip_items + length(iter_terms)
  }

  # Interaction terms
  intervars <- if (is.null(interactionvars)) X else interactionvars
  if (!is.null(notinteractionvars)) {
    intervars <- intervars[!intervars %in% notinteractionvars]
  }
  inter.dt <- combn(intervars, 2) |> t() |> data.table()
  inters <- inter.dt[, paste(V1, V2, sep = ":")]

  terms <- inters
  if (!is.null(quadraticvars)) terms <- c(quad_terms, terms)
  if (quadratic_interactions) terms <- c(terms, quad_inters)
  #if (!is.null(itervars)) terms <- c(iter_terms, terms)

  inter_terms <- NULL
  for (i in seq_along(terms)) {
    if (i %in% skip_items) next
    if (is.null(inter_terms)) {
      inter_terms <- terms[i]
    } else {
      idx <- length(inter_terms)
      inter_terms <- c(
        inter_terms,
        paste(inter_terms[length(inter_terms)], terms[i], sep = "+")
      )
    }
  }

  # random effects
  random_effects <- match.arg(random_effects)
  if (is.null(id)) id <- "ID"
  if (is.null(slopevars)) slopevars <- "T"
  re <- if(random_effects == "intercepts") {
    sprintf("(1|%s)", id)
  } else if (random_effects == "slopes") {
    sprintf("(%s|%s)", paste(slopevars, collapse = "+"), id)
  } else {
    NULL
  }

  formulas <- if (is.null(itervars)) {
    X <- paste(X, collapse = "+")
    c(sprintf("%s ~ %s", Y, X), sprintf("%s ~ %s+%s", Y, X, inter_terms))
  } else {
    Xb <- X[!X %in% itervars]
    Xb <- paste(Xb, collapse = "+")
    X <- paste(X, collapse = "+")
    c(sprintf("%s ~ %s", Y, Xb),
      sprintf("%s ~ %s+%s", Y, Xb, iter_terms),
      sprintf("%s ~ %s+%s", Y, X, inter_terms))
  }

  if (!is.null(re)) formulas <- paste(formulas, re, sep = "+")
  #if (!is.null(skip_items)) {
    #for (i in seq_along(skip_items))
      #j <- skip_items[i] + 1 - i
      #formulas <- formulas[-j]
  #}

  return(formulas)
}
