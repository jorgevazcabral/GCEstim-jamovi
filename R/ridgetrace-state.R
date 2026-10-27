# Retire only obsolete presentation options; keep the version-1 numerical cache.
# jmvcore sends NONE tombstones to the host so these options are removed on save.
.gceRTCleanOptions <- function(options) {
  for (key in c("plotHeatmap", "plotSelected", "normalize"))
    if (options$has(key)) options$.removeOption(key)
}

.gceRTCleanImages <- function(results) {
  allowed <- c("axisScale", "showIntercept", "extrema", "facet", "displayTerms", "seBars", "foldMean")
  for (key in c("plotPaths", "plotRanges", "plotSensitivity", "plotCV", "plotSample", "plotFolds", "plotStability")) {
    image <- results[[key]]
    state <- image$state
    if (is.list(state$options)) {
      cleaned <- state$options[intersect(names(state$options), allowed)]
      if (!identical(cleaned, state$options)) {
        state$options <- cleaned
        image$setState(state)
      }
    }
  }
}

# Plain numerical state; estimation is exclusively GCEstim::ridgetrace().
.gceRTLambda <- function(custom, text, minimum, maximum, count) {
  if (!isTRUE(custom)) {
    if (any(!is.finite(c(minimum, maximum, count))) || minimum <= 0 ||
        maximum <= minimum || count < 2 || count != floor(count))
      jmvcore::reject("Use positive finite lambda limits with minimum < maximum and at least two grid values.", code = "invalid_lambda_grid")
    return(NULL)
  }
  parts <- trimws(strsplit(paste0(text, ",END"), ",", fixed = TRUE)[[1L]])
  values <- suppressWarnings(as.numeric(head(parts, -1L)))
  if (length(values) < 2L || any(!is.finite(values)) || any(values <= 0) || anyDuplicated(values))
    jmvcore::reject("Custom lambda values must contain at least two distinct, strictly positive finite numbers separated by commas.", code = "invalid_lambda_values")
  sort(values)
}

.gceRTCompact <- function(fit, p, cv, folds) {
  lambda <- as.numeric(fit$lambda); b <- fit$coef.lambda
  if (length(lambda) < 2L || any(!is.finite(lambda)) || any(lambda <= 0) ||
      is.unsorted(lambda, strictly = TRUE) || !is.matrix(b) ||
      !identical(dim(b), c(as.integer(p), length(lambda))) ||
      is.null(rownames(b)) || length(fit$error.lambda) != length(lambda) ||
      any(vapply(fit[c("min.coef", "max.coef", "max.abs.coef")], length, integer(1)) != p))
    stop("GCEstim returned incompatible ridge path dimensions or lambda values.")
  if (any(!is.finite(b))) stop("GCEstim returned nonfinite coefficient estimates.")
  if (isTRUE(cv)) {
    e <- as.matrix(fit$error.lambda.cv); bc <- fit$coef.lambda.cv
    if (!identical(dim(e), c(as.integer(folds), length(lambda))) || length(bc) != folds ||
        !all(vapply(bc, function(x) is.matrix(x) && identical(dim(x), dim(b)) &&
          identical(rownames(x), rownames(b)), logical(1))))
      stop("GCEstim returned incompatible cross-validation dimensions.")
  } else { e <- NULL; bc <- NULL }
  # Rounded matrix column names are never used for indexing.
  colnames(b) <- NULL
  if (!is.null(e)) colnames(e) <- NULL
  bc <- lapply(bc, function(x) { colnames(x) <- NULL; x })
  list(lambda = lambda, coefficients = b,
       minimum = as.numeric(fit$min.coef), maximum = as.numeric(fit$max.coef),
       maxAbs = as.numeric(fit$max.abs.coef), maxResidual = as.numeric(fit$max.abs.residual),
       sample = as.numeric(fit$error.lambda), cvErrors = e,
       foldCoefficients = if (isTRUE(cv)) bc else NULL)
}

.gceRTSelection <- function(fit, rule, manual, cv) {
  L <- length(fit$lambda); mean <- se <- rep(NA_real_, L)
  if (isTRUE(cv)) {
    valid <- apply(is.finite(fit$cvErrors), 2L, all)
    if (any(valid)) {
      mean[valid] <- colMeans(fit$cvErrors[, valid, drop = FALSE])
      se[valid] <- apply(fit$cvErrors[, valid, drop = FALSE], 2L, stats::sd) / sqrt(nrow(fit$cvErrors))
    }
    errors <- mean
  } else errors <- fit$sample
  valid <- which(is.finite(errors))
  minimum <- if (length(valid)) valid[which.min(errors[valid])] else NA_integer_
  oneSE <- NA_integer_
  if (isTRUE(cv) && !is.na(minimum) && is.finite(se[minimum])) {
    candidates <- which(is.finite(mean) & mean <= mean[minimum] + se[minimum])
    oneSE <- candidates[which.max(fit$lambda[candidates])]
  }
  warning <- character()
  if (any(!is.finite(fit$sample)) || (isTRUE(cv) && any(!is.finite(fit$cvErrors))))
    warning <- "Some returned errors are undefined or nonfinite. CV mean/SE requires all folds to be finite; unavailable errors are excluded from automatic selection."
  if (rule == "1se" && !isTRUE(cv))
    jmvcore::reject("The 1-SE rule requires cross-validation. Choose Minimum error (in-sample) or Manual evaluated lambda.", code = "one_se_requires_cv")
  selected <- if (rule == "1se") oneSE else minimum
  if (rule == "manual") {
    if (!is.finite(manual) || manual <= 0)
      jmvcore::reject("Select a positive finite lambda from the evaluated grid.", code = "invalid_manual_lambda")
    matches <- which(fit$lambda == manual)
    if (length(matches) != 1L)
      jmvcore::reject("Manual lambda must match exactly one evaluated grid value. No interpolation or refitting is performed.", code = "manual_lambda_not_evaluated")
    selected <- matches
  }
  if (is.na(selected)) warning <- c(warning,
    "No finite error-based lambda selection is available. Coefficient paths and grid summaries remain available; use Manual evaluated lambda to inspect coefficients.")
  list(index = selected, minimum = minimum, oneSE = oneSE, mean = mean, se = se, warning = warning)
}

.gceRTUpsert <- function(table, key, values) {
  if (key %in% unlist(table$rowKeys)) table$setRow(rowKey = key, values = values)
  else table$addRow(rowKey = key, values = values)
}

.gceRTTerms <- function(fit, showIntercept, text) {
  terms <- rownames(fit$coefficients)
  if (!isTRUE(showIntercept)) terms <- setdiff(terms, "(Intercept)")
  if (nzchar(trimws(text))) {
    wanted <- trimws(scan(text = text, what = character(), sep = ",", quote = '"', quiet = TRUE))
    unknown <- setdiff(wanted[nzchar(wanted)], rownames(fit$coefficients))
    if (length(unknown)) stop(paste("Unknown coefficient names:", paste(unknown, collapse = ", ")))
    terms <- intersect(terms, wanted)
  }
  if (!length(terms)) stop("No coefficients are selected for this plot.")
  terms
}

.gceRTLong <- function(fit, terms) {
  b <- fit$coefficients[terms, , drop = FALSE]
  data.frame(term = factor(rep(terms, length(fit$lambda)), levels = terms),
             lambda = rep(fit$lambda, each = length(terms)), estimate = as.numeric(b))
}

# Presentation maps are keyed by the complete internal coefficient universe.
.gceRTColours <- function(fit, theme) {
  terms <- unique(rownames(fit$coefficients))
  terms <- terms[order(sub("^`([^`]*)`", "\\1", terms), terms, method = "radix")]
  # theme$color contains two foreground colours, not the coefficient palette.
  palette <- if (!is.null(theme$palette))
    jmvcore::colorPalette(6L, theme$palette, "color") else theme$color
  if (is.null(palette) || !length(palette))
    palette <- jmvcore::colorPalette(6L, "jmv", "color")
  # Cycle the existing theme palette deterministically for additional coefficients.
  stats::setNames(rep_len(palette, length(terms)), terms)
}

.gceRTLabels <- function(terms) {
  terms <- unique(as.character(terms))
  labels <- sub("^`([^`]*)`", "\\1", terms)
  # Disambiguate globally, including names that already resemble a suffix.
  labels <- make.unique(labels, sep = " [coefficient] ")
  stats::setNames(labels, terms)
}

.gceRTLineType <- function(rule) {
  switch(rule %||% "min", min = "dotted", `1se` = "dashed", manual = "solid",
    stop("Unknown ridge selection rule."))
}

.gceRTPlot <- function(cache, state, ggtheme, theme) {
  f <- cache$fit; s <- state$selection; key <- state$key; o <- state$options; lambda <- f$lambda
  coefficientPlot <- key %in% c("plotPaths", "plotRanges", "plotSensitivity", "plotStability")
  if (coefficientPlot) {
    labels <- .gceRTLabels(rownames(f$coefficients))
    terms <- .gceRTTerms(f, o$showIntercept, o$displayTerms); d <- .gceRTLong(f, terms)
  }
  if (!is.null(f$cvErrors)) {
    summary <- .gceRTSelection(f, "min", lambda[1L], TRUE)
    s$mean <- summary$mean; s$se <- summary$se
  }
  lambdaAxis <- TRUE
  selectedLambda <- if (!is.na(s$index)) lambda[s$index] else NA_real_
  colour <- (theme$color %||% "#346da4")[[1L]]
  if (key == "plotPaths") {
    p <- ggplot2::ggplot(d, ggplot2::aes(lambda, estimate, colour = term, group = term)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey65", linetype = "dashed") +
      ggplot2::geom_line() + ggplot2::labs(y = "Coefficient (original units)", colour = NULL)
    if (isTRUE(o$extrema)) {
      extrema <- function(maximum) do.call(rbind, lapply(terms, function(term) {
        x <- d[d$term == term, , drop = FALSE]
        x[if (maximum) which.max(x$estimate) else which.min(x$estimate), , drop = FALSE]
      }))
      p <- p + ggplot2::geom_point(data = extrema(TRUE),
        ggplot2::aes(x = lambda, y = estimate, colour = term),
        inherit.aes = FALSE, shape = 24, size = 1) +
        ggplot2::geom_point(data = extrema(FALSE),
          ggplot2::aes(x = lambda, y = estimate, colour = term),
          inherit.aes = FALSE, shape = 25, size = 1)
    }
    if (isTRUE(o$facet)) p <- p + ggplot2::facet_wrap(~term, scales = "free_y", labeller = ggplot2::as_labeller(labels))
  } else if (key == "plotCV") {
    d <- data.frame(lambda = lambda, error = s$mean, se = s$se)
    p <- ggplot2::ggplot(d, ggplot2::aes(lambda, error)) +
      ggplot2::geom_line(colour = colour, na.rm = TRUE) + ggplot2::geom_point(colour = colour, na.rm = TRUE) +
      ggplot2::labs(y = paste("Mean CV", cache$measure),
        caption = if (isTRUE(o$seBars)) "Bars: mean +/- 1 SE across folds; not confidence intervals for generalization error." else NULL)
    if (isTRUE(o$seBars)) p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = error - se, ymax = error + se), width = 0, na.rm = TRUE)
  } else if (key == "plotSample") {
    d <- data.frame(lambda = lambda, error = replace(f$sample, !is.finite(f$sample), NA_real_))
    p <- ggplot2::ggplot(d, ggplot2::aes(lambda, error)) +
      ggplot2::geom_line(colour = colour, na.rm = TRUE) + ggplot2::geom_point(colour = colour, na.rm = TRUE) +
      ggplot2::labs(y = paste("In-sample", cache$measure), caption = "In-sample fit does not estimate optimal predictive performance.")
  } else if (key == "plotRanges") {
    positions <- match(terms, rownames(f$coefficients))
    d <- data.frame(term = factor(terms, levels = rev(terms)),
      minimum = f$minimum[positions], maximum = f$maximum[positions],
      estimate = if (is.na(s$index)) NA_real_ else f$coefficients[positions, s$index])
    p <- ggplot2::ggplot(d, ggplot2::aes(x = estimate, y = term, colour = term)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey65", linetype = "dashed")
    p <- p + ggplot2::geom_segment(
      ggplot2::aes(x = minimum, xend = maximum, yend = term))
    p <- p + ggplot2::geom_point(size = 2, na.rm = TRUE) +
      ggplot2::labs(x = "Coefficient (original units)", y = NULL, colour = NULL,
        caption = paste0("Ranges over evaluated grid; not confidence intervals. Points: selected coefficients. ",
          "Selected lambda: ", if (is.finite(selectedLambda)) format(selectedLambda, digits = 8) else "unavailable")) +
      ggplot2::guides(colour = "none")
    lambdaAxis <- FALSE
  } else if (key == "plotSensitivity") {
    b <- f$coefficients[terms, , drop = FALSE]
    sensitivity <- (b[, -1L, drop = FALSE] - b[, -ncol(b), drop = FALSE]) /
      rep(diff(log10(lambda)), each = length(terms))
    sensitivity[!is.finite(sensitivity)] <- NA_real_
    d <- data.frame(term = factor(rep(terms, length(lambda)-1L), levels = terms),
      lambda = rep(10^(head(log10(lambda), -1) + diff(log10(lambda))/2), each = length(terms)),
      estimate = as.numeric(sensitivity))
    p <- ggplot2::ggplot(d, ggplot2::aes(lambda, estimate, colour = term, group = term)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65") +
      ggplot2::geom_line(na.rm = TRUE) + ggplot2::geom_point(na.rm = TRUE) +
      ggplot2::facet_wrap(~term, scales = "free_y", labeller = ggplot2::as_labeller(labels)) +
      ggplot2::labs(y = "Coefficient change / change in log10(lambda)", colour = NULL,
        caption = "Finite differences at log-interval midpoints: descriptive and grid-dependent; not significance or an optimal-lambda criterion.")
  } else if (key == "plotFolds") {
    e <- f$cvErrors
    if (is.null(e)) stop("Cross-validation is disabled.")
    e[!is.finite(e)] <- NA_real_
    d <- data.frame(lambda = rep(lambda, each = nrow(e)), error = as.numeric(e),
      fold = factor(rep(seq_len(nrow(e)), length(lambda))))
    p <- ggplot2::ggplot(d, ggplot2::aes(lambda, error, group = fold)) +
      ggplot2::geom_line(colour = colour, alpha = .35, na.rm = TRUE) +
      ggplot2::labs(y = paste("Fold CV", cache$measure), caption = "Thin lines: held-out fold errors.")
    if (isTRUE(o$foldMean)) p <- p + ggplot2::geom_line(data = data.frame(lambda = lambda, error = s$mean),
      ggplot2::aes(lambda, error), inherit.aes = FALSE, linewidth = .9, colour = colour, na.rm = TRUE)
  } else if (key == "plotStability") {
    if (is.null(f$foldCoefficients)) stop("Cross-validation is disabled.")
    folds <- do.call(rbind, lapply(seq_along(f$foldCoefficients), function(k) {
      values <- f$foldCoefficients[[k]][terms, , drop = FALSE]
      values[!is.finite(values)] <- NA_real_
      data.frame(term = factor(rep(terms, length(lambda)), levels = terms),
        lambda = rep(lambda, each = length(terms)), estimate = as.numeric(values), fold = k)
    }))
    p <- ggplot2::ggplot(d, ggplot2::aes(lambda, estimate, colour = term, group = term)) +
      ggplot2::geom_line(data = folds, ggplot2::aes(group = interaction(term, fold), linetype = "Fold training"),
        alpha = .35, na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(linetype = "Full data"), linewidth = .9) +
      ggplot2::facet_wrap(~term, scales = "free_y", labeller = ggplot2::as_labeller(labels)) +
      ggplot2::scale_linetype_manual(values = c("Fold training" = "dotted", "Full data" = "solid")) +
      ggplot2::labs(y = "Coefficient (original units)", colour = NULL, linetype = NULL,
        caption = "Fold-training variation; not bootstrap distributions or coefficient confidence intervals.")
  } else stop("Unknown ridge plot.")
  if (key == "plotRanges") p <- p + ggplot2::scale_y_discrete(labels = labels)
  if (lambdaAxis && key != "plotCV" && is.finite(selectedLambda)) p <- p + ggplot2::geom_vline(xintercept = selectedLambda, linetype = .gceRTLineType(s$rule), colour = "grey40")
  if (key == "plotCV") {
    marks <- data.frame(lambda = lambda[c(s$minimum, s$oneSE)], rule = c("Minimum error", "1-SE rule"))
    marks <- marks[is.finite(marks$lambda), , drop = FALSE]
    if (nrow(marks)) p <- p + ggplot2::geom_vline(data = marks,
      ggplot2::aes(xintercept = lambda, linetype = rule), colour = "black") +
      ggplot2::scale_linetype_manual(values = c("Minimum error" = .gceRTLineType("min"), "1-SE rule" = .gceRTLineType("1se")),
        breaks = c("Minimum error", "1-SE rule")) + ggplot2::labs(linetype = NULL)
  }
  if (lambdaAxis && o$axisScale == "log") p <- p + ggplot2::scale_x_log10()
  if (lambdaAxis) p <- p + ggplot2::labs(x = if (o$axisScale == "log") "Lambda (log10 scale)" else "Lambda")
  if (!is.null(ggtheme)) p <- p + ggtheme
  # jamovi themes include colour scales; install the stable named scale last.
  if (coefficientPlot) p <- p + ggplot2::scale_colour_manual(values = .gceRTColours(f, theme), labels = labels)
  if (!is.null(p$labels$caption)) p <- p + ggplot2::labs(caption = paste(strwrap(p$labels$caption, width = 80), collapse = "\n"))
  if (coefficientPlot && (key %in% c("plotSensitivity", "plotStability") || (key == "plotPaths" && isTRUE(o$facet))))
    p <- p + ggplot2::guides(colour = "none")
  else if (key == "plotPaths") p <- p + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 3))
  p + ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 8), legend.position = "bottom", legend.box = "vertical")
}

.gceRTRender <- function(image, cache, ggtheme, theme) {
  if (is.null(image$state) || is.null(cache$fit)) return(FALSE)
  tryCatch({ print(.gceRTPlot(cache, image$state, ggtheme, theme)); TRUE }, error = function(e) FALSE)
}
