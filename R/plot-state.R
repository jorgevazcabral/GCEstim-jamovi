# Plot states contain plain data; fitted models and ggplot environments stay out of images.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Validate optional plot input without rejecting the fitted analysis.
.gceTrueCoefficients <- function(text, form, data) {
  tryCatch({
    text <- trimws(text)
    if (!nzchar(text))
      stop("The precision-error plots require a vector of true coefficients.")
    parts <- trimws(strsplit(paste0(text, ",END"), ",", fixed = TRUE)[[1L]])
    values <- suppressWarnings(as.numeric(head(parts, -1L)))
    if (!length(values) || any(!is.finite(values)))
      stop("True coefficients must be finite numeric values separated by commas.")
    expected <- ncol(stats::model.matrix(form, data))
    if (length(values) != expected)
      stop(sprintf("The true coefficient vector must contain %d values for the current model (including the intercept when selected).", expected))
    list(values = values, warning = NULL)
  }, error = function(e) list(values = NULL,
    warning = paste0("Precision-error plots skipped: ", conditionMessage(e))))
}

# Only coefficient vectors and plotting metadata are retained, never a model.
# Reconstitute the small interface used by GCEstim's precision plot method on
# demand; changing true coefficients then requires no estimation or CV.
.gcePrecisionBasis <- function(best) {
  fields <- c("coefficients", "error", "support.ok", "support.stdUL",
              "support.signal.min", "support.signal.1se", "support.signal.elbow",
              "support.signal.manual", "nep")
  basis <- lapply(fields, function(key) best[[key]])
  names(basis) <- fields
  sr <- best$results$nocv$support.results
  candidates <- lapply(sr[as.character(best$support.ok)], function(x)
    list(coefficients = as.numeric(x$coefficients)))
  # These three summary entries are expected by plot.lmgce(which = 5).
  if (length(sr) == 1L)
    candidates <- list(list(coefficients = as.numeric(sr[[1L]]$coefficients)))
  for (key in if (length(sr) > 1L) tail(names(sr), 3L) else character())
    candidates[[key]] <- sr[[key]]
  basis$results <- list(
    nocv = list(support.results = candidates),
    twosteps = lapply(best$results$twosteps, function(x)
      list(coefficients = as.numeric(x$coefficients))),
    OLS = list(error = as.numeric(best$results$OLS$error),
               matrix.coef = best$results$OLS$matrix.coef))
  basis
}

.gcePrecisionState <- function(basis, which, trueCoef) {
  p <- plot(structure(basis, class = "lmgce"), which = which,
            coef = trueCoef, OLS = TRUE, NormEnt = TRUE)
  if (is.list(p) && !inherits(p, "ggplot")) p <- p[[1L]]
  .gceDiagnosticState(p, which,
    supportLines = if (which == 5L) .gceSupportLines(basis) else NULL,
    selectedSupport = if (which == 5L) basis$support.stdUL else NULL,
    selectedNEP = if (which == 5L) basis$nep else NULL)
}

# Image states hold a reference and presentation settings, not a second copy of
# cached data. Legacy self-contained image states can still render on reopening.
.gcePlotImage <- function(image, cache) {
  state <- image$state
  if (is.null(state$key)) {
    if (identical(state$which, 4L))
      return(list(state = .gceEntropyState(state, cache$plots$plot3)))
    return(image)
  }
  presentation <- state
  state <- cache$plots[[presentation$key]]
  if (is.null(state)) return(list(state = NULL))
  if (presentation$key %in% c("plot1", "plot3", "plot4"))
    state$coefficientTerms <- rownames(cache$coefficients) %||% cache$plots$plot1$data$term
  if (presentation$key == "plot4")
    state <- .gceEntropyState(state, cache$plots$plot3)
  if (!isTRUE(presentation$OLS)) {
    if (presentation$key == "plot1") state$data$ols <- NA_real_
    if (presentation$key == "plot3") state$ols <- NULL
    state$overlays <- Filter(function(x) !"yintercept" %in% names(x), state$overlays)
  }
  if (presentation$key %in% c("plot2", "plot5"))
    state$showNEP <- isTRUE(presentation$NormEnt)
  list(state = state)
}

# Recovery diagnostics belong in results, never in exported graphical images.
.gceEntropyWarning <- function(results, options, states) {
  state <- if (isTRUE(options$plot4)) .gceEntropyState(states$plot4, states$plot3) else NULL
  message <- state$entropyNotice
  notice <- jmvcore::Notice$new(options = options, name = ".entropyWarning", type = jmvcore::NoticeType$WARNING)
  notice$setContent(message %||% "")
  notice$setVisible(!is.null(message) && nzchar(message))
  results$insert(1, notice)
}

.gceCompactCVFit <- function(fit) {
  best <- fit$best
  list(
    best = list(
      support.matrix = as.matrix(best$support.matrix),
      p = as.matrix(best$p),
      w = as.matrix(best$w),
      fitted.values = as.numeric(best$fitted.values),
      residuals = as.numeric(best$residuals)
    ),
    results = as.data.frame(fit$results),
    support.signal.points.best = fit$support.signal.points.best
  )
}

# Generate all available plot data while the selected GCEstim model exists.
# The returned states contain data frames and atomic values only.
.gceCachePlotStates <- function(fit, ci, level, method,
                                supportGridAvailable,
                                reestimationPlotsAvailable,
                                cvPlotAvailable) {
  best <- fit$best
  warnings <- character()
  safely <- function(label, expr) tryCatch(expr, error = function(e) {
    warnings <<- c(warnings, paste0(label, ": ", conditionMessage(e)))
    NULL
  })
  states <- list(
    plot1 = safely("Estimates and confidence intervals",
                   .gceEstimateState(best, ci, level, method, includeOLS = TRUE)),
    plotCV = if (cvPlotAvailable) safely("Cross-validation error",
                                        .gceCVState(fit)) else NULL
  )
  ols <- safely("OLS comparison", .gceOLSState(best))
  for (which in c(2L, 3L, 4L, 6L)) {
    if ((which %in% 2L:5L && !supportGridAvailable) ||
        (which %in% 6L:7L && !reestimationPlotsAvailable))
      next
    key <- paste0("plot", which)
    states[[key]] <- safely(key, {
      plotFit <- best

      p <- plot(plotFit, which = which, coef = NULL,
                OLS = TRUE, NormEnt = TRUE)
      if (is.list(p) && !inherits(p, "ggplot")) p <- p[[1L]]
      .gceDiagnosticState(
        p, which, showNEP = FALSE,
        supportLines = if (which %in% 2L:5L) .gceSupportLines(best) else NULL,
        ols = if (which == 3L) ols else NULL,
        selectedSupport = if (which %in% c(2L, 5L)) best$support.stdUL else NULL,
        selectedError = if (which == 2L) best$error.measure.cv.mean else NULL,
        selectedNEP = if (which == 2L) best$nep.cv.mean else if (which == 5L)
          best$nep else NULL
      )
    })
  }
  if (!is.null(states$plot4)) {
    states$plot4$coefficientTerms <- names(stats::coef(best))
    states$plot4 <- .gceEntropyState(states$plot4, states$plot3)
  }
  list(states = states, warnings = warnings)
}


# plot3 retains full-data entropy rows that the upstream entropy plot can drop
# because of NA in an unused ratio. Recover presentation only from that same
# cached stage; never estimate or calculate entropy here.
.gceEntropyState <- function(state, source = NULL) {
  if (is.null(state)) return(NULL)
  fields <- c("support", "predictor", "nepk")
  clean <- function(data) {
    if (!is.data.frame(data) || !all(fields %in% names(data)) ||
        !is.numeric(data$support) || !is.numeric(data$nepk)) return(NULL)
    ids <- as.character(data$predictor)
    data$predictor <- ids
    data <- data[is.finite(data$support) & data$support > 0 &
                   !is.na(ids) & nzchar(ids) & is.finite(data$nepk), , drop = FALSE]
    if (anyDuplicated(data[c("support", "predictor")])) return(NULL)
    data
  }
  available <- clean(state$data)
  candidate <- if (identical(source$which, 3L)) clean(source$data) else NULL
  terms <- state$coefficientTerms %||% unique(as.character(state$data$predictor))
  compatible <- !is.null(candidate) && nrow(candidate) > 0 &&
    all(candidate$predictor %in% terms) &&
    (is.null(source$supportLines) || is.null(state$supportLines) ||
       identical(source$supportLines, state$supportLines))
  if (compatible && !is.null(available) && nrow(available)) {
    overlap <- merge(available[fields], candidate[fields],
                     by = c("support", "predictor"), suffixes = c(".old", ".new"))
    compatible <- nrow(overlap) == nrow(available) &&
      all(abs(overlap$nepk.old - overlap$nepk.new) <= 1e-12)
  }
  if (compatible) {
    state$data <- candidate
    state$entropyComplete <- TRUE
    state$entropyNotice <- NULL
  } else {
    state$data <- available %||% data.frame(support = numeric(), predictor = character(), nepk = numeric())
    if (!isTRUE(state$entropyComplete))
      state$entropyNotice <- paste(
        "Complete entropy trajectories are unavailable in this saved state.",
        "Only available observations are shown; missing values have not been estimated.")
  }
  state
}

# Wrap coefficient legends within the existing image width.
.gceCoefficientLegend <- function(terms) {
  labels <- .gceCoefficientLabels(terms)
  if (!length(labels)) return(list(labels = labels, columns = 1L))
  columns <- max(1L, min(length(labels),
    max(ceiling(length(labels) / 5), floor(58 / (max(nchar(labels)) + 7)))))
  width <- max(8L, floor(58 / columns) - 7L)
  wrap <- function(label) {
    words <- strsplit(label, " ", fixed = TRUE)[[1L]]
    words <- unlist(lapply(words, function(word) {
      starts <- seq.int(1L, max(1L, nchar(word)), by = width)
      substring(word, starts, starts + width - 1L)
    }), use.names = FALSE)
    paste(strwrap(paste(words, collapse = " "), width = width), collapse = "\n")
  }
  list(labels = stats::setNames(vapply(labels, wrap, character(1)), names(labels)),
       columns = as.integer(columns),
       textSize = if (length(labels) > 5L || max(nchar(labels)) > 28L) 8 else NULL)
}

.gceEntropyLegend <- function(terms) .gceCoefficientLegend(terms)

# Complete internal identifiers define presentation; never use displayed labels as keys.
.gceCoefficientLabels <- function(terms) {
  terms <- unique(as.character(terms))
  labels <- sub("^`([^`]*)`", "\\1", terms)
  stats::setNames(make.unique(labels, sep = " [coefficient] "), terms)
}

.gceTermColours <- function(terms, theme) {
  terms <- unique(as.character(terms))
  # theme$color is the foreground pair, not jamovi's categorical palette.
  colours <- if (!is.null(theme$palette))
    jmvcore::colorPalette(max(6L, length(terms)), theme$palette, "color") else theme$color
  if (is.null(colours) || !length(colours))
    colours <- jmvcore::colorPalette(max(6L, length(terms)), "jmv", "color")
  stats::setNames(rep_len(colours, length(terms)), terms)
}

.gceCoefficientTerms <- function(state, terms) {
  state$coefficientTerms %||% unique(as.character(terms))
}

.gceLabelTables <- function(results, terms, supportTerms = terms, probabilityTerms = terms) {
  labels <- .gceCoefficientLabels(unique(c(terms, supportTerms, probabilityTerms)))
  for (key in c("coefficients", "signalSupports", "signalProbabilities")) {
    ids <- switch(key, coefficients = terms, signalSupports = supportTerms,
                  signalProbabilities = probabilityTerms)
    table <- results[[key]]; keys <- unlist(table$rowKeys)
    if (length(keys) != length(ids)) next
    for (i in seq_along(ids))
      table$setRow(rowKey = keys[i], values = list(term = unname(labels[ids[i]])))
  }
}

# GCEstim's selected-point frames retain coefficient identifiers as row names.
# Read those identities for legacy states too; never guess by sorted position.
.gceCoefficientOverlay <- function(overlay, terms) {
  ids <- if ("predictor" %in% names(overlay)) as.character(overlay$predictor) else rownames(overlay)
  if (length(ids) == nrow(overlay) && !anyDuplicated(ids) && all(ids %in% terms))
    overlay$predictor <- ids
  overlay
}

.gceEstimateState <- function(fit, ci, level, method, includeOLS = FALSE) {
  estimates <- stats::coef(fit)
  terms <- names(estimates)
  if (is.null(terms))
    terms <- rownames(ci)
  
  intervals <- as.matrix(ci)[terms, , drop = FALSE]
  
  data <- data.frame(
    term = terms,
    estimate = as.numeric(estimates),
    lower = as.numeric(intervals[, 1L]),
    upper = as.numeric(intervals[, 2L]),
    stringsAsFactors = FALSE
  )
  
  data$ols <- NA_real_
  
  if (isTRUE(includeOLS)) {
    olsFit <- fit$results$OLS$res
    
    if (!is.null(olsFit)) {
      olsCoef <- stats::coef(olsFit)
      
      if (all(terms %in% names(olsCoef)))
        data$ols <- as.numeric(olsCoef[terms])
    }
  }
  
  list(
    data = data,
    coefficientTerms = terms,
    caption = sprintf(
      "Estimate (%s%% CI, %s)", 100 * level,
      switch(method, z = "asymptotic", percentile = "bootstrap percentile",
             basic = "bootstrap basic", method)
    )
  )
}

.gceRenderEstimates <- function(image, ggtheme, theme) {
  state <- image$state
  if (is.null(state) || is.null(state$data) || !nrow(state$data))
    return(FALSE)
  data <- state$data
  # Older saved image states used an internal CI code in their caption.
  # Translate it when rendering so reopening an .omv needs no model refit.
  caption <- state$caption %||% 'Estimate'
  methods <- c(z = "asymptotic", percentile = "bootstrap percentile",
               basic = "bootstrap basic")
  for (method in names(methods))
    caption <- sub(paste0("^Estimates \\(", method, " ([0-9.]+)% CI\\)$"),
                   paste0("Estimate (\\1% CI, ", methods[[method]], ")"), caption)
  terms <- .gceCoefficientTerms(state, data$term)
  labels <- .gceCoefficientLabels(terms)
  data$term <- factor(data$term, levels = data$term)
  p <- ggplot2::ggplot(data, ggplot2::aes(x = estimate, y = term)) +
    ggplot2::geom_vline(xintercept = 0,
                        linetype = 'dashed',
                        colour = 'grey50') +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = lower,
        xmax = upper,
        colour = term
      ),
      orientation = 'y',
      width = 0
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = term), size = 2) +
    ggplot2::labs(
      x = caption,
      y = NULL,
      colour = NULL
    )
  if (any(is.finite(data$ols))) {
    p <- p + ggplot2::geom_point(
      data = data[is.finite(data$ols), , drop = FALSE],
      ggplot2::aes(x = .data$ols, y = .data$term, colour = .data$term),
      inherit.aes = FALSE,
      shape = 1,
      size = 2
    )
  }
  if (!is.null(ggtheme))
    p <- p + ggtheme
  p <- p + ggplot2::scale_colour_manual(values = .gceTermColours(terms, theme), labels = labels) +
    ggplot2::scale_y_discrete(labels = labels)
  p <- p + ggplot2::guides(colour = "none") +
    ggplot2::theme(legend.position = "none")
  print(p)
  TRUE
}

.gceCVState <- function(fit) {
  results <- as.data.frame(fit$results)
  list(
    data = data.frame(
      signalPoints = results$support.signal.points,
      noisePoints = results$support.noise.points,
      weight = results$weight,
      cvMean = results$error.measure.cv.mean
    ),
    measure = as.character(fit$best$error)
  )
}

.gceRenderCV <- function(image, ggtheme, theme) {
  state <- image$state
  if (is.null(state) || is.null(state$data) || !nrow(state$data))
    return(FALSE)
  data <- state$data
  data$noisePoints <- factor(data$noisePoints)
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = signalPoints,
      y = cvMean,
      group = noisePoints,
      colour = noisePoints
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap( ~ weight, ncol = 1, scales = 'free') +
    ggplot2::labs(x = 'Number of points of the signal support',
                  y = paste0('CV-', state$measure),
                  colour = NULL)
  if (!is.null(theme$color) && length(theme$color))
    p <- p + ggplot2::scale_colour_manual(values = rep_len(theme$color, nlevels(data$noisePoints)))
  if (!is.null(ggtheme))
    p <- p + ggtheme
  p <- p + ggplot2::theme(legend.position = 'bottom',
                          legend.title = ggplot2::element_blank())
  print(p)
  TRUE
}

.gceRidgeState <- function(ridgeFit) {
  # GCEstim's plot method supplies the long-format trace and marked extrema.
  p <- plot(ridgeFit, log = TRUE, range = TRUE)
  extrema <- lapply(p$layers, function(layer)
    layer$data)
  extrema <- Filter(function(x)
    is.data.frame(x) &&
      all(
        c('variable', 'max_x', 'max_y', 'min_x', 'min_y') %in% names(x)
      ), extrema)
  list(data = as.data.frame(p$data),
       coefficientTerms = rownames(ridgeFit$coef.lambda),
       extrema = if (length(extrema))
         extrema[[1L]]
       else
         NULL)
}

.gceRenderRidge <- function(image, ggtheme, theme) {
  state <- image$state
  if (is.null(state) || is.null(state$data) || !nrow(state$data))
    return(FALSE)
  p <- ggplot2::ggplot(
    state$data,
    ggplot2::aes(
      x = lambda,
      y = coefficient,
      group = variable,
      colour = variable
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = 'Log ridge penalty',
      y = 'Coefficients',
      colour = NULL
    )
  if (!is.null(state$extrema)) {
    p <- p + ggplot2::geom_point(
      data = state$extrema,
      ggplot2::aes(x = max_x, y = max_y, colour = variable),
      inherit.aes = FALSE,
      shape = 24,
      size = 1
    ) +
      ggplot2::geom_point(
        data = state$extrema,
        ggplot2::aes(x = min_x, y = min_y, colour = variable),
        inherit.aes = FALSE,
        shape = 25,
        size = 1
      )
  }
  if (!is.null(ggtheme))
    p <- p + ggtheme
  terms <- .gceCoefficientTerms(state, state$data$variable)
  legend <- .gceCoefficientLegend(terms)
  p <- p + ggplot2::scale_colour_manual(values = .gceTermColours(terms, theme),
    labels = legend$labels) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = legend$columns, byrow = TRUE))
  if (!is.null(legend$textSize))
    p <- p + ggplot2::theme(legend.text = ggplot2::element_text(size = legend$textSize))
  p <- p + ggplot2::theme(legend.position = 'bottom',
                          legend.title = ggplot2::element_blank())
  print(p)
  TRUE
}

.gceSupportLines <- function(fit) {
  list(
    min = fit$support.signal.min,
    oneSE = fit$support.signal.1se,
    manual = fit$support.signal.manual
  )
}

.gceOLSState <- function(fit) {
  olsFit <- fit$results$OLS$res
  if (is.null(olsFit))
    return(NULL)
  
  values <- stats::coef(olsFit)
  terms <- names(stats::coef(fit))
  
  if (!all(terms %in% names(values)))
    return(NULL)
  
  data.frame(
    predictor = terms,
    estimate = as.numeric(values[terms])
  )
}

# GCEstim computes the diagnostic statistics. Keep their plotted rows and the
# small reference-layer coordinates, then construct the graph on jamovi's device.
.gceDiagnosticState <- function(plot,
                                which,
                                showNEP = FALSE,
                                supportLines = NULL,
                                ols = NULL,
                                selectedSupport = NULL,
                                selectedError = NULL,
                                selectedNEP = NULL) {
  if (is.null(plot) || !inherits(plot, 'ggplot'))
    return(NULL)
  overlays <- lapply(plot$layers, function(layer) {
    data <- layer$data
    if (!is.data.frame(data) || !nrow(data))
      return(NULL)
    keep <- intersect(c('x', 'y', 'yintercept'), names(data))
    if (!length(keep))
      return(NULL)
    overlay <- as.data.frame(data[, unique(c(keep, intersect("predictor", names(data)))), drop = FALSE])
    if (which %in% c(3L, 4L))
      overlay <- .gceCoefficientOverlay(overlay, unique(as.character(plot$data$predictor)))
    overlay
  })
  overlays <- Filter(Negate(is.null), overlays)
  
  if (which == 5L && is.null(selectedError)) {
    for (layer in plot$layers) {
      layerColour <- layer$aes_params$colour %||% layer$aes_params$color
      layerData <- layer$data
      
      if (is.character(layerColour) && length(layerColour) == 1L &&
          grepl("red", layerColour, ignore.case = TRUE) &&
          is.data.frame(layerData) && nrow(layerData) == 1L &&
          all(c("x", "y") %in% names(layerData)) &&
          isTRUE(all.equal(
            as.numeric(layerData$x),
            as.numeric(selectedSupport),
            tolerance = 1e-7
          ))) {
        selectedError <- as.numeric(layerData$y)
        break
      }
    }
  }
  
  label <- function(x)
    if (is.character(x) && length(x) == 1L)
      x
  else
    NULL
  # Keep only the columns read by .gceRenderDiagnostic. In particular, do not
  # retain arbitrary attributes or model objects attached to a plot's data.
  plotData <- as.data.frame(plot$data)
  plotData <- plotData[, intersect(
    c("support", "predictor", "estimate", "nepk", "error", "nep",
      "error.LL", "error.UL", "nep.LL", "nep.UL"), names(plotData)
  ), drop = FALSE]
  list(
    data = plotData,
    coefficientTerms = if (which %in% c(3L, 4L)) unique(as.character(plot$data$predictor)) else NULL,
    which = as.integer(which),
    overlays = overlays,
    showNEP = isTRUE(showNEP),
    supportLines = supportLines,
    ols = ols,
    selectedSupport = selectedSupport,
    selectedError = selectedError,
    selectedNEP = selectedNEP,
    xLabel = label(plot$labels$x),
    yLabel = label(plot$labels$y),
    title = label(plot$labels$title)
  )
}

.gceRenderDiagnostic <- function(image, ggtheme, theme) {
  state <- image$state
  if (identical(state$which, 4L)) state <- .gceEntropyState(state)
  if (is.null(state) || is.null(state$data) ||
      (!nrow(state$data) && !identical(state$which, 4L)))
    return(FALSE)
  
  data <- state$data
  which <- state$which
  
  colour <- if (!is.null(theme$color) && length(theme$color))
    theme$color[[1L]]
  else
    "#346da4"
  
  entropyColour <- if (!is.null(theme$color) &&
                       length(theme$color) > 1L)
    theme$color[[2L]]
  else
    "forestgreen"
  
  hasEntropyAxis <- FALSE
  
  if (which %in% c(3L, 4L)) {
    y <- if (which == 3L)
      "estimate"
    else
      "nepk"
    if (!all(c("support", "predictor", y) %in% names(data)))
      return(FALSE)
    
    data$value <- data[[y]]
    data$predictor <- factor(data$predictor, levels = sort(unique(as.character(data$predictor))))
    
    p <- ggplot2::ggplot(data,
                         ggplot2::aes(
                           x = support,
                           y = value,
                           group = predictor,
                           colour = predictor
                         )) +
      ggplot2::geom_line() +
      ggplot2::labs(
        x = state$xLabel %||% "Flexibility parameter",
        y = state$yLabel %||% if (which == 3L)
          "Estimates"
        else
          "Normalized Entropy",
        colour = NULL
      )
    
    if (which == 3L && !is.null(state$ols)) {
      ols <- state$ols
      ols <- ols[is.finite(ols$estimate) &
                   ols$predictor %in% levels(data$predictor), , drop = FALSE]
      
      if (nrow(ols))
        p <- p + ggplot2::geom_hline(
          data = ols,
          ggplot2::aes(yintercept = estimate, colour = predictor),
          inherit.aes = FALSE,
          linetype = "dotted"
        )
    }
    
  } else {
    if (!all(c("support", "error") %in% names(data)))
      return(FALSE)
    
    supportPlot <- which %in% c(2L, 5L)
    
    selectedPoint <- function(value) {
      x <- state$selectedSupport
      
      if (!supportPlot ||
          !is.numeric(x) || length(x) != 1L || !is.finite(x) ||
          !is.numeric(value) || length(value) != 1L ||
          !is.finite(value))
        return(NULL)
      
      data.frame(support = x, value = value)
    }
    
    p <- ggplot2::ggplot(data, ggplot2::aes(x = support, y = error)) +
      ggplot2::geom_line(colour = colour,
                         linetype = if (supportPlot)
                           "solid"
                         else
                           "dashed") +
      ggplot2::labs(
        x = state$xLabel %||% if (which %in% c(6L, 7L))
          "Number of re-estimations"
        else
          "Flexibility parameter",
        y = state$yLabel %||% if (which %in% c(2L, 6L))
          "Prediction error"
        else
          "Precision error"
      )
    
    errorPoint <- selectedPoint(state$selectedError)
    if (!is.null(errorPoint))
      p <- p + ggplot2::geom_point(
        data = errorPoint,
        ggplot2::aes(y = value),
        colour = "black",
        size = 2
      )
    
    if (!supportPlot)
      p <- p + ggplot2::geom_point(colour = colour, size = 2)
    
    if (all(c("error.LL", "error.UL") %in% names(data)))
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = error.LL, ymax = error.UL),
        width = 0,
        colour = "grey50"
      )
    
    if (isTRUE(state$showNEP) &&
        all(c("nep", "error") %in% names(data))) {
      left <- range(data$error, finite = TRUE)
      right <- range(data$nep, finite = TRUE)
      
      if (all(is.finite(c(left, right))) &&
          diff(right) > 0 && diff(left) > 0) {
        slope <- diff(left) / diff(right)
        shift <- left[[1L]] - slope * right[[1L]]
        
        data$nepScaled <- shift + slope * data$nep
        
        if (all(c("nep.LL", "nep.UL") %in% names(data))) {
          data$nepLowerScaled <- shift + slope * data$nep.LL
          data$nepUpperScaled <- shift + slope * data$nep.UL
        }
        
        p <- p +
          ggplot2::geom_line(
            data = data,
            ggplot2::aes(y = nepScaled),
            colour = entropyColour,
            linetype = "solid"
          ) +
          ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ (. - shift) / slope, name = "Normalized Entropy"))
        
        if (all(c("nepLowerScaled", "nepUpperScaled") %in% names(data)))
          p <- p + ggplot2::geom_errorbar(
            data = data,
            ggplot2::aes(
              ymin = nepLowerScaled,
              ymax = nepUpperScaled
            ),
            width = 0,
            colour = entropyColour
          )
        
        hasEntropyAxis <- TRUE
        
        nepPoint <- selectedPoint(state$selectedNEP)
        if (!is.null(nepPoint)) {
          nepPoint$value <- shift + slope * nepPoint$value
          
          p <- p + ggplot2::geom_point(
            data = nepPoint,
            ggplot2::aes(y = value),
            colour = entropyColour,
            size = 2
          )
        }
      }
    }
    
    if (which %in% c(6L, 7L))
      p <- p + ggplot2::scale_x_continuous(breaks = data$support)
  }
  
  for (overlay in state$overlays) {
    if (all(c("x", "y") %in% names(overlay)) &&
        !which %in% c(2L, 5L)) {
      if (which %in% c(3L, 4L))
        overlay <- .gceCoefficientOverlay(overlay, .gceCoefficientTerms(state, data$predictor))
      if (which %in% c(3L, 4L) && "predictor" %in% names(overlay) &&
          all(overlay$predictor %in% .gceCoefficientTerms(state, data$predictor))) {
        
        p <- p + ggplot2::geom_point(
          data = overlay,
          ggplot2::aes(
            x = x,
            y = y,
            colour = predictor
          ),
          inherit.aes = FALSE,
          size = 2
        )
        
      } else {
        p <- p + ggplot2::geom_point(
          data = overlay,
          ggplot2::aes(x = x, y = y),
          inherit.aes = FALSE,
          colour = colour,
          size = 2
        )
      }
      
    } else if ("yintercept" %in% names(overlay) &&
               which != 3L) {
      p <- p + ggplot2::geom_hline(
        yintercept = as.numeric(overlay$yintercept),
        linetype = "dotted",
        colour = "grey55"
      )
    }
  }
  
  if (which %in% 2L:5L && length(state$supportLines)) {
    lineTypes <- c(min = "dotted",
                   oneSE = "dashed",
                   manual = "dotted")
    lineColours <- c(
      min = "black",
      oneSE = "black",
      manual = "green4"
    )
    
    for (name in intersect(names(lineTypes), names(state$supportLines))) {
      value <- state$supportLines[[name]]
      
      if (is.numeric(value) &&
          length(value) == 1L &&
          is.finite(value))
        p <- p + ggplot2::geom_vline(xintercept = value,
                                     linetype = lineTypes[[name]],
                                     colour = lineColours[[name]])
    }
  }
  
  if (!is.null(ggtheme))
    p <- p + ggtheme
  
  if (which %in% c(3L, 4L)) {
    terms <- .gceCoefficientTerms(state, unique(as.character(state$data$predictor)))
    legend <- .gceCoefficientLegend(terms)
    p <- p + ggplot2::scale_colour_manual(values = .gceTermColours(terms, theme),
      labels = legend$labels,
      limits = if (length(terms)) NULL else character())
    {
      p <- p + ggplot2::guides(colour = ggplot2::guide_legend(ncol = legend$columns, byrow = TRUE))
      if (!is.null(legend$textSize))
        p <- p + ggplot2::theme(legend.text = ggplot2::element_text(size = legend$textSize))
    }
  }
  if (which %in% 2L:5L)
    p <- p + ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_blank()
    )
  
  if (hasEntropyAxis)
    p <- p + ggplot2::theme(
      axis.title.y.right = ggplot2::element_text(colour = entropyColour),
      axis.text.y.right = ggplot2::element_text(colour = entropyColour),
      axis.ticks.y.right = ggplot2::element_line(colour = entropyColour),
      axis.line.y.right = ggplot2::element_line(colour = entropyColour)
    )
  
  print(p)
  TRUE
}
