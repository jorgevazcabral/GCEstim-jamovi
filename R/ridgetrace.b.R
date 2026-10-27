RidgeTraceClass <- R6::R6Class(
  "RidgeTraceClass", inherit = RidgeTraceBase,
  private = list(
    .run = function() {
      .gceRTCleanOptions(self$options)
      .gceRTCleanImages(self$results)
      if (!isTRUE(self$options$run)) {
        self$results$text$setVisible(TRUE)
        self$results$text$setContent("<p>Select variables and options, then select <b>Fit and update model</b>.</p>")
        return()
      }
      o <- self$options
      dep <- unlist(o$dep); covs <- unlist(o$covs); factors <- unlist(o$factors)
      if (length(dep) != 1L || !length(c(covs, factors)))
        jmvcore::reject("Select one dependent variable and at least one covariate or factor.", code = "variables_required")
      if (anyDuplicated(c(dep, covs, factors)))
        jmvcore::reject("Use distinct variables for the dependent variable, covariates and factors.", code = "duplicate_variables")
      data <- self$data[, c(dep, covs, factors), drop = FALSE]
      data <- data[stats::complete.cases(data), , drop = FALSE]
      if (nrow(data) < 3L)
        jmvcore::reject("This analysis requires at least 3 complete cases.", code = "too_few_cases")
      if (!all(vapply(data[, c(dep, covs), drop = FALSE], function(x) is.numeric(x) && all(is.finite(x)), logical(1))))
        jmvcore::reject("The dependent variable and covariates must contain finite numeric values.", code = "nonfinite_data")
      if (length(unique(data[[dep]])) < 2L)
        jmvcore::reject("The dependent variable has no variation.", code = "constant_response")
      for (name in factors) data[[name]] <- droplevels(as.factor(data[[name]]))
      constant <- c(covs[vapply(data[, covs, drop = FALSE], function(x) length(unique(x)) < 2L, logical(1))],
        factors[vapply(data[, factors, drop = FALSE], function(x) nlevels(x) < 2L, logical(1))])
      if (length(constant)) jmvcore::reject("The following predictors have no variation: {variables}.",
        code = "constant_predictors", variables = paste(constant, collapse = ", "))
      formStr <- jmvcore::constructFormula(dep = dep, terms = c(covs, factors))
      if (!isTRUE(o$intercept)) formStr <- paste(formStr, "- 1")
      form <- stats::as.formula(formStr)
      X <- tryCatch(stats::model.matrix(form, data), error = function(e)
        jmvcore::reject("Could not construct the model: {message}", code = "invalid_model", message = conditionMessage(e)))
      if (isTRUE(o$standardize) && !isTRUE(o$intercept))
        jmvcore::reject("Standardization without an intercept is not currently supported by the GCEstim implementation. This is an implementation limitation, not a statistical restriction.", code = "unsupported_standardized_no_intercept")
      if (isTRUE(o$standardize) && isTRUE(o$intercept) && isTRUE(o$cv) && ncol(X) == 2L)
        jmvcore::reject("Cross-validation with standardization, an intercept and one predictor column is not currently supported by the GCEstim implementation. This is an implementation limitation, not a statistical restriction.", code = "unsupported_single_predictor_cv")
      if (isTRUE(o$cv) && (o$cvNfolds < 2L || o$cvNfolds > nrow(data)))
        jmvcore::reject("The number of CV folds must be between 2 and the number of complete cases.", code = "invalid_folds")
      lambda <- .gceRTLambda(o$customGrid, o$lambdaValues, o$lambdaMin, o$lambdaMax, o$lambdaN)
      cached <- self$results$fitCache$state
      reused <- is.list(cached) && identical(cached$version, 1L) && is.list(cached$fit)
      self$results$text$setVisible(TRUE)
      self$results$text$setContent(if (reused) "<p>Updating cached ridge results...</p>" else "<p>Computing ridge trace...</p>")
      private$.checkpoint()
      if (!reused) {
        started <- Sys.time(); warnings <- character()
        fit <- tryCatch(withCallingHandlers(GCEstim::ridgetrace(
          formula = form, data = data, lambda = lambda,
          lambda.min = o$lambdaMin, lambda.max = o$lambdaMax, lambda.n = o$lambdaN,
          standardize = o$standardize, penalize.intercept = o$penalizeIntercept,
          cv = o$cv, cv.nfolds = o$cvNfolds, seed = o$seed, errormeasure = o$errorMeasure),
          warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
          error = function(e) jmvcore::reject("Ridge estimation failed in GCEstim: {message}", code = "ridge_fit_failed", message = conditionMessage(e)))
        compact <- tryCatch(.gceRTCompact(fit, ncol(X), o$cv, o$cvNfolds),
          error = function(e) jmvcore::reject("Invalid ridge results: {message}", code = "invalid_ridge_results", message = conditionMessage(e)))
        cached <- list(version = 1L, fit = compact, measure = o$errorMeasure,
          engineVersion = as.character(utils::packageVersion("GCEstim")), warnings = warnings,
          elapsed = as.numeric(difftime(Sys.time(), started, units = "secs")))
        self$results$fitCache$setState(cached)
      }
      f <- cached$fit
      if (o$selection == "manual" && (!is.finite(o$manualIndex) || o$manualIndex < 1L ||
          o$manualIndex != floor(o$manualIndex) || o$manualIndex > length(f$lambda)))
        jmvcore::reject("Manual lambda grid position must identify an evaluated lambda (from 1 to {count}).",
          code = "invalid_manual_index", count = length(f$lambda))
      manual <- if (o$selection == "manual") f$lambda[o$manualIndex] else f$lambda[1L]
      s <- .gceRTSelection(f, o$selection, manual, o$cv)
      selected <- if (is.na(s$index)) rep(NA_real_, nrow(f$coefficients)) else f$coefficients[, s$index]
      selectedLambda <- if (is.na(s$index)) NA_real_ else f$lambda[s$index]
      rule <- if (o$selection == "manual") "Manual evaluated lambda" else if (o$selection == "1se") "1-SE rule" else if (o$cv) "Minimum mean CV error" else "Minimum in-sample error (descriptive)"
      summary <- c(formula = paste(deparse(form), collapse = " "), n = nrow(data),
        predictors = ncol(X) - as.integer(o$intercept), intercept = if (o$intercept) "Included" else "Omitted",
        standardize = if (o$standardize) "Yes" else "No",
        penalize = if (!o$intercept) "Not applicable" else if (o$standardize) "No (intercept reconstructed after standardization)" else if (o$penalizeIntercept) "Yes" else "No",
        grid = if (o$customGrid) "Custom, sorted ascending" else "Geometric",
        minimum = min(f$lambda), maximum = max(f$lambda), count = length(f$lambda),
        cv = if (o$cv) "Enabled" else "Disabled", folds = if (o$cv) o$cvNfolds else "Not applicable",
        seed = if (o$cv) o$seed else "Not applicable", error = cached$measure,
        rule = rule, selected = if (is.finite(selectedLambda)) format(selectedLambda, digits = 15) else "Unavailable", engine = cached$engineVersion)
      labels <- c("Formula", "Complete cases", "Predictor columns (excluding intercept)", "Intercept", "Standardization", "Intercept penalization",
        "Lambda grid", "Minimum lambda", "Maximum lambda", "Lambda values", "Cross-validation", "CV folds", "Seed", "Error measure", "Selection rule", "Selected lambda", "GCEstim version")
      for (i in seq_along(summary)) .gceRTUpsert(self$results$modelSummary, names(summary)[i], list(measure = labels[i], value = as.character(summary[i])))
      self$results$coefficientRanges$setNote(key = "descriptive", note = "Ranges describe the evaluated lambda grid; they are not confidence intervals.", init = TRUE)
      coefficientLabels <- .gceRTLabels(rownames(f$coefficients))
      for (j in seq_len(nrow(f$coefficients))) {
        term <- rownames(f$coefficients)[j]
        .gceRTUpsert(self$results$coefficients, term, list(term = unname(coefficientLabels[term]), estimate = as.numeric(selected[j])))
        .gceRTUpsert(self$results$coefficientRanges, term, list(term = unname(coefficientLabels[term]),
          minimum = f$minimum[j], maximum = f$maximum[j], maxAbs = f$maxAbs[j], selected = as.numeric(selected[j])))
      }
      self$results$coefficients$setTitle(paste("Coefficients at Lambda", if (is.finite(selectedLambda)) format(selectedLambda, digits = 8) else "Unavailable"))
      self$results$lambdaResults$setVisible(isTRUE(o$showLambdaTable))
      self$results$lambdaResults$getColumn("sample")$setTitle(paste("In-sample", cached$measure))
      for (column in c("cvMean", "cvSE", "oneSE")) self$results$lambdaResults$getColumn(column)$setVisible(isTRUE(o$cv))
      if (isTRUE(o$showLambdaTable)) {
        for (i in seq_len(min(length(f$lambda), o$tableRows))) .gceRTUpsert(self$results$lambdaResults, as.character(i), list(
          position = i, lambda = f$lambda[i], sample = if (is.finite(f$sample[i])) f$sample[i] else NA_real_,
          cvMean = s$mean[i], cvSE = s$se[i], selected = if (!is.na(s$index) && i == s$index) "Yes" else "",
          minimum = if (!is.na(s$minimum) && i == s$minimum) "Yes" else "",
          oneSE = if (!is.na(s$oneSE) && i == s$oneSE) "Yes" else ""))
        self$results$lambdaResults$setNote(key = "limit", note = paste("Showing up to", o$tableRows, "rows in ascending lambda order. CV summaries require all folds to be finite. SE describes fold variation, not a confidence interval."), init = TRUE)
      }
      plotWarnings <- character()
      optionNames <- c("axisScale", "showIntercept", "extrema", "facet", "displayTerms", "seBars", "foldMean")
      plotOptions <- stats::setNames(lapply(optionNames, function(key) o[[key]]), optionNames)
      for (key in c("plotPaths", "plotCV", "plotSample", "plotRanges", "plotSensitivity", "plotFolds", "plotStability")) {
        available <- if (key %in% c("plotCV", "plotFolds", "plotStability")) isTRUE(o$cv) else TRUE
        visible <- available && if (key == "plotSample") {
          o$plotSample == "show" || (o$plotSample == "auto" && !isTRUE(o$cv))
        } else isTRUE(o[[key]])
        image <- self$results[[key]]; image$setVisible(visible)
        if (!visible) next
        state <- list(presentationVersion = 1L, key = key, selection = c(s[c("index", "minimum", "oneSE")], list(rule = o$selection)), options = plotOptions)
        ok <- tryCatch({ ggplot2::ggplot_build(.gceRTPlot(cached, state, NULL, list())); TRUE }, error = function(e) {
          plotWarnings <<- c(plotWarnings, paste0(key, " skipped: ", conditionMessage(e))); FALSE })
        if (key %in% c("plotPaths", "plotRanges", "plotSensitivity", "plotStability") && ok) {
          terms <- .gceRTTerms(f, o$showIntercept, o$displayTerms)
          height <- if (key == "plotRanges") max(400, 120 + 24 * length(terms)) else if (key == "plotPaths" && !isTRUE(o$facet)) max(400, 320 + 24 * ceiling(length(terms)/3)) else max(400, 200 * ceiling(length(terms)/3))
          image$setSize(650, height)
        }
        image$setState(if (ok) state else NULL)
        if (!ok) image$setVisible(FALSE)
      }
      notes <- c(if (isTRUE(o$cv) && isTRUE(o$standardize))
        "GCEstim standardizes using the complete dataset before constructing CV folds. These CV results preserve full-data preprocessing and do not use independent fold-specific preprocessing.",
        if (!isTRUE(o$cv)) "Minimum in-sample error describes fit to these observations; it does not estimate optimal predictive performance.",
        if (isTRUE(o$cv)) "Fold-based SE values are descriptive, not confidence intervals for generalization error.")
      notice <- jmvcore::Notice$new(options = o, name = ".method", type = jmvcore::NoticeType$INFO)
      notice$setContent(paste(notes, collapse = " ")); self$results$insert(1, notice)
      warnings <- unique(c(cached$warnings, s$warning, plotWarnings,
        if (length(f$foldCoefficients) && any(vapply(f$foldCoefficients, function(x) any(!is.finite(x)), logical(1))))
          "Some fold-training coefficients are nonfinite and are omitted from stability plots."))
      if (length(warnings)) {
        notice <- jmvcore::Notice$new(options = o, name = ".fitWarning", type = jmvcore::NoticeType$WARNING)
        notice$setContent(paste(warnings, collapse = "; ")); self$results$insert(1, notice)
      }
      notice <- jmvcore::Notice$new(options = o, name = ".fitTime", type = jmvcore::NoticeType$INFO)
      notice$setContent(paste0(if (reused) "Previous fit reused; original computation time: " else "Computation time: ", round(cached$elapsed, 3), " s."))
      self$results$insert(1, notice); self$results$text$setVisible(FALSE)
    },
    .plotRanges = function(image, ggtheme, theme) .gceRTRender(image, self$results$fitCache$state, ggtheme, theme),
    .plotSensitivity = function(image, ggtheme, theme) .gceRTRender(image, self$results$fitCache$state, ggtheme, theme),
    .plotFolds = function(image, ggtheme, theme) .gceRTRender(image, self$results$fitCache$state, ggtheme, theme),
    .plotStability = function(image, ggtheme, theme) .gceRTRender(image, self$results$fitCache$state, ggtheme, theme),
    .plotPaths = function(image, ggtheme, theme) .gceRTRender(image, self$results$fitCache$state, ggtheme, theme),
    .plotCV = function(image, ggtheme, theme) .gceRTRender(image, self$results$fitCache$state, ggtheme, theme),
    .plotSample = function(image, ggtheme, theme) .gceRTRender(image, self$results$fitCache$state, ggtheme, theme)
  )
)
