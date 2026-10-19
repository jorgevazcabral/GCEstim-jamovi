# Plot states contain plain data; fitted models and ggplot environments stay out of images.
`%||%` <- function(x, y) if (is.null(x)) y else x

.gcePackFitCache <- function(fit, ci, coefficients, plots) {
  bytes <- serialize(list(version = 3L, fit = fit, ci = ci,
                          coefficients = coefficients, plots = plots), NULL)
  as.integer(memCompress(bytes, type = "gzip"))
}

.gceUnpackFitCache <- function(payload) {
  if (is.null(payload))
    return(NULL)
  
  tryCatch(
    unserialize(memDecompress(
      as.raw(as.integer(unlist(payload, use.names = FALSE))),
      type = "gzip"
    )),
    error = function(e) NULL
  )
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
.gceCachePlotStates <- function(fit, ci, level, method, trueCoef,
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
  for (which in 2L:7L) {
    if ((which %in% 2L:5L && !supportGridAvailable) ||
        (which %in% 6L:7L && !reestimationPlotsAvailable) ||
        (which %in% c(5L, 7L) && is.null(trueCoef)))
      next
    key <- paste0("plot", which)
    states[[key]] <- safely(key, {
      plotFit <- best
      if (which == 5L) {
        # GCEstim's precision plot uses every support fit, but its x values
        # include only successful supports. Align them if a candidate failed.
        sr <- best$results$nocv$support.results
        plotFit$results$nocv$support.results <-
          sr[c(as.character(best$support.ok), tail(names(sr), 3L))]
      }
      p <- plot(plotFit, which = which, coef = trueCoef,
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
  list(states = states, warnings = warnings)
}

.gceTermColours <- function(terms, theme) {
  terms <- sort(unique(as.character(terms)))
  colours <- theme$color
  if (is.null(colours) || !length(colours))
    colours <- grDevices::hcl.colors(length(terms), 'Dark 3')
  stats::setNames(rep_len(colours, length(terms)), terms)
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
  p <- p + ggplot2::scale_colour_manual(values = .gceTermColours(data$term, theme))
  if (!is.null(ggtheme))
    p <- p + ggtheme
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
  p <- p + ggplot2::scale_colour_manual(values = .gceTermColours(state$data$variable, theme))
  if (!is.null(ggtheme))
    p <- p + ggtheme
  p <- p + ggplot2::guides(colour = "none") +
    ggplot2::theme(legend.position = "none")
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
    as.data.frame(data[, keep, drop = FALSE])
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
  if (is.null(state) || is.null(state$data) || !nrow(state$data))
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
      ) +
      ggplot2::scale_colour_manual(values = .gceTermColours(data$predictor, theme))
    
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
      if (which %in% c(3L, 4L) &&
          nrow(overlay) == nlevels(data$predictor)) {
        overlay$predictor <- levels(data$predictor)
        
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
