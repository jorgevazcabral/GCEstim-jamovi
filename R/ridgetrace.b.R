RidgeTraceClass <- R6::R6Class(
  "RidgeTraceClass", inherit = RidgeTraceBase,
  private = list(
    .run = function() {
      if (!isTRUE(self$options$run)) {
        self$results$text$setVisible(TRUE)
        self$results$text$setContent("<p>Select variables and options, then select <b>Fit and update model</b>.</p>")
        return()
      }
      o <- self$options
      effective <- .gceRTEffective(o)
      settings <- .gceRTSettings(o, effective)
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
      # Reject invalid presentation input once, before fitting or constructing plots.
      previous <- self$results$fitCache$state
      displayFit <- if (is.list(previous) && identical(previous$version, 1L) && is.list(previous$fit) &&
        (is.null(previous$settings) || identical(previous$settings, settings))) previous$fit else list(coefficients = t(X))
      displayTerms <- .gceRTTerms(displayFit, o$showIntercept, o$displayTerms)
      if (effective$standardize && isTRUE(o$cv) && ncol(X) == 2L)
        jmvcore::reject("Cross-validation with standardization, an intercept and one predictor column is not currently supported by the GCEstim implementation. This is an implementation limitation, not a statistical restriction.", code = "unsupported_single_predictor_cv")
      if (isTRUE(o$cv) && (o$cvNfolds < 2L || o$cvNfolds > nrow(data)))
        jmvcore::reject("The number of CV folds must be between 2 and the number of complete cases.", code = "invalid_folds")
      lambda <- .gceRTLambda(o$customGrid, o$lambdaValues, o$lambdaMin, o$lambdaMax, o$lambdaN)
      cached <- self$results$fitCache$state
      reused <- is.list(cached) && identical(cached$version, 1L) && is.list(cached$fit)
      if (reused && !is.null(cached$settings)) reused <- identical(cached$settings, settings)
      self$results$text$setVisible(TRUE)
      self$results$text$setContent(if (reused) "<p>Updating cached ridge results...</p>" else "<p>Computing ridge trace...</p>")
      private$.checkpoint()
      if (!reused) {
        started <- Sys.time(); warnings <- character()
        fit <- tryCatch(withCallingHandlers(GCEstim::ridgetrace(
          formula = form, data = data, lambda = lambda,
          lambda.min = o$lambdaMin, lambda.max = o$lambdaMax, lambda.n = o$lambdaN,
          standardize = effective$standardize, penalize.intercept = o$penalizeIntercept,
          cv = o$cv, cv.nfolds = o$cvNfolds, seed = o$seed, errormeasure = o$errorMeasure),
          warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
          error = function(e) jmvcore::reject("Ridge estimation failed in GCEstim: {message}", code = "ridge_fit_failed", message = conditionMessage(e)))
        compact <- tryCatch(.gceRTCompact(fit, ncol(X), o$cv, o$cvNfolds),
          error = function(e) jmvcore::reject("Invalid ridge results: {message}", code = "invalid_ridge_results", message = conditionMessage(e)))
        cached <- list(version = 1L, fit = compact, measure = o$errorMeasure,
          engineVersion = as.character(utils::packageVersion("GCEstim")), warnings = warnings,
          elapsed = as.numeric(difftime(Sys.time(), started, units = "secs")), settings = settings)
        self$results$fitCache$setState(cached)
      }
      f <- cached$fit
      if (o$selection == "manual" && (!is.finite(o$manualIndex) || o$manualIndex < 1L ||
          o$manualIndex != floor(o$manualIndex) || o$manualIndex > length(f$lambda)))
        jmvcore::reject("Manual lambda grid position must identify an evaluated lambda (from 1 to {count}).",
          code = "invalid_manual_index", count = length(f$lambda))
      manual <- if (o$selection == "manual") f$lambda[o$manualIndex] else f$lambda[1L]
      s <- .gceRTSelection(f, effective$selection, manual, o$cv)
      selected <- if (is.na(s$index)) rep(NA_real_, nrow(f$coefficients)) else f$coefficients[, s$index]
      selectedLambda <- if (is.na(s$index)) NA_real_ else f$lambda[s$index]
      rule <- if (effective$selection == "manual") "Manual (grid position)" else if (effective$selection == "1se") "1-SE rule" else if (o$cv) "Minimum mean CV error" else "Minimum in-sample error (descriptive)"
      summary <- c(formula = paste(deparse(form), collapse = " "), n = nrow(data),
        predictors = ncol(X) - as.integer(o$intercept), intercept = if (o$intercept) "Included" else "Omitted",
        standardize = if (!o$intercept) "Not applicable" else if (effective$standardize) "Yes" else "No",
        penalize = if (!o$intercept) "Not applicable" else if (effective$standardize) "No (intercept reconstructed after standardization)" else if (o$penalizeIntercept) "Yes" else "No",
        grid = if (o$customGrid) "Custom, sorted ascending" else "Geometric",
        count = length(f$lambda),
        cv = if (o$cv) "Enabled" else "Disabled", folds = if (o$cv) o$cvNfolds else "Not applicable",
        seed = if (o$cv) o$seed else "Not applicable", error = cached$measure,
        rule = rule, engine = cached$engineVersion)
      labels <- c("Formula", "Complete cases", "Predictor columns (excluding intercept)", "Intercept", "Standardization", "Intercept penalization",
        "Lambda grid", "Lambda values", "Cross-validation", "CV folds", "Seed", "Error measure", "Selection rule", "GCEstim version")
      for (i in seq_along(summary)) .gceRTUpsert(self$results$modelSummary, names(summary)[i], list(measure = labels[i], value = as.character(summary[i])))
      .gceRTUpsert(self$results$selectedLambda, "selection", list(minimum = min(f$lambda), maximum = max(f$lambda), selected = selectedLambda, rule = rule))
      self$results$coefficientRanges$setNote(key = "descriptive", note = "Ranges describe the evaluated lambda grid; they are not confidence intervals.", init = TRUE)
      coefficientLabels <- .gceCoefficientLabels(rownames(f$coefficients))
      for (j in seq_len(nrow(f$coefficients))) {
        term <- rownames(f$coefficients)[j]
        .gceRTUpsert(self$results$coefficients, term, list(term = unname(coefficientLabels[term]), estimate = as.numeric(selected[j])))
        .gceRTUpsert(self$results$coefficientRanges, term, list(term = unname(coefficientLabels[term]),
          minimum = f$minimum[j], maximum = f$maximum[j], maxAbs = f$maxAbs[j], selected = as.numeric(selected[j])))
      }
      self$results$coefficients$setTitle(paste("Coefficients at Lambda", if (is.finite(selectedLambda)) format(signif(selectedLambda, 3), digits = 3, trim = TRUE) else "Unavailable"))
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
      for (key in c("plotPaths", "plotCV", "plotSample", "plotRanges", "plotSensitivity", "plotFolds", "plotStability")) {
        available <- if (key %in% c("plotCV", "plotFolds", "plotStability")) isTRUE(o$cv) else TRUE
        visible <- available && if (key == "plotSample") {
          o$plotSample == "show" || (o$plotSample == "auto" && !isTRUE(o$cv))
        } else isTRUE(o[[key]])
        image <- self$results[[key]]; image$setVisible(visible)
        if (!visible) next
        optionNames <- .gceRTPlotOptions(key)
        plotOptions <- stats::setNames(lapply(optionNames, function(name) o[[name]]), optionNames)
        selection <- if (key == "plotCV") s[c("minimum", "oneSE")] else if (key == "plotRanges") s["index"] else c(s["index"], list(rule = effective$selection))
        state <- list(presentationVersion = 2L, key = key, selection = selection, options = plotOptions)
        if (key %in% c("plotPaths", "plotRanges", "plotSensitivity", "plotStability")) state$terms <- displayTerms
        if (reused && identical(image$state, state) && image$status != "error") next
        image$setStatus("complete")
        ok <- tryCatch({ ggplot2::ggplot_build(.gceRTPlot(cached, state, NULL, list())); TRUE }, error = function(e) {
          state$error <<- paste0(image$title, ": ", conditionMessage(e))
          image$setError(state$error); FALSE })
        if (key %in% c("plotPaths", "plotRanges", "plotSensitivity", "plotStability") && ok) {
          terms <- displayTerms
          height <- if (key == "plotRanges") max(400, 120 + 24 * length(terms)) else if (key == "plotPaths" && !isTRUE(o$facet)) max(400, 320 + 24 * ceiling(length(terms)/3)) else max(400, 200 * ceiling(length(terms)/3))
          image$setSize(650, height)
        }
        image$setState(state)
      }
      notes <- c(if (isTRUE(o$cv) && effective$standardize)
        "GCEstim standardizes using the complete dataset before constructing CV folds. These CV results preserve full-data preprocessing and do not use independent fold-specific preprocessing.",
        if (!isTRUE(o$cv) && o$selection == "1se") "CV is disabled; minimum in-sample error is being used for descriptive selection. The stored 1-SE preference is retained.")
      notice <- jmvcore::Notice$new(options = o, name = ".method", type = jmvcore::NoticeType$INFO)
      notice$setContent(paste(notes, collapse = " ")); notice$setVisible(length(notes) > 0L); self$results$insert(1, notice)
      warnings <- unique(c(cached$warnings, s$warning,
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
