TASWGCEClass <- R6::R6Class(
  "TASWGCEClass",
  inherit = TASWGCEBase,
  private = list(
    
    .init = function() { 
      coefTable <- self$results$coefficients
      level <- format(100 * self$options$bootConfLevel, trim = TRUE)
      
      method <- if (isTRUE(self$options$bootstrap)) {
        paste("bootstrap", self$options$bootCIMethod)
      } else {
        "asymptotic"
      }
      
      ciTitle <- paste0(level, "% Confidence Interval (", method, ")")
      coefTable$getColumn("ciLower")$setSuperTitle(ciTitle)
      coefTable$getColumn("ciUpper")$setSuperTitle(ciTitle)
      
      if (isTRUE(self$options$bootstrap)) {
        coefTable$setNote(
          key = "bootstrapColumns",
          note = "Asymptotic SE, z, and p are hidden when bootstrap intervals are selected.",
          init = TRUE
        )
      }
      
      summaryRows <- c(
        formula = "Formula",
        n = "Complete cases",
        cvNfolds = "Cross-validation folds",
        errorMeasure = "Prediction-error measure",
        errorMeasureWhich = "Selection rule",
        seed = "Seed",
        method = "Optimization method",
        supportMethod = "Support method",
        signalSupportSpecification = "Signal support specification",
        noiseSupportMethod = "Noise support specification",
        supportSignalVectorN = "Number of support spaces",
        M = "Candidate signal support points",
        J = "Candidate noise support points",
        weight = "Candidate noise weight values",
        twostepsN = "Post-GCE re-estimations",
        ciType = "Confidence interval type",
        ciLevel = "Confidence level"
      )
      
      if (isTRUE(self$options$bootstrap)) {
        summaryRows <- c(
          summaryRows,
          bootB = "Bootstrap samples",
          bootMethod = "Bootstrap resampling method"
        )
      }
      
      weightValues <- suppressWarnings(as.numeric(
        trimws(strsplit(self$options$weight, ",", fixed = TRUE)[[1L]])
      ))
      weightValues <- unique(weightValues[is.finite(weightValues)])
      
      for (plotName in paste0("plot", 1:5)) {
        width <- switch(plotName,
                        plot1 = 520,
                        plot2 = 560,
                        plot3 = 520,
                        plot4 = 520,
                        plot5 = 560
        )
        self$results[[plotName]]$setSize(width, 400)
      }
      
      self$results$plotCV$setSize(
        700,
        160 + 140 * max(1L, length(weightValues))
      )

      summaryTable <- self$results$modelSummary
      for (key in names(summaryRows)) {
        if (!(key %in% unlist(summaryTable$rowKeys))) {
          summaryTable$addRow(
            rowKey = key,
            values = list(measure = summaryRows[[key]])
          )
        }
      }
      
    },
    
    .run = function() {
        
      dep  <- unlist(self$options$dep)
      covs    <- unlist(self$options$covs)
      factors <- unlist(self$options$factors)

      for (plotName in paste0("plot", 1:5)) {
        width <- switch(plotName,
                        plot1 = 520,
                        plot2 = 560,
                        plot3 = 520,
                        plot4 = 520,
                        plot5 = 560
        )
        self$results[[plotName]]$setSize(width, 400)
      }
      
      self$results$text$setVisible(TRUE)
      
      if (!isTRUE(self$options$run)) {
        self$results$text$setContent(
          "<p>Select variables and options, then select <b>Fit and update model</b>.</p>"
        )
        return()
      }
      
      if (length(dep) == 0 || (length(covs) == 0 && length(factors) == 0)) {
        jmvcore::reject(
          "Select one dependent variable and at least one covariate or factor.",
          code = "variables_required"
        )
      }
      
      supportGridAvailable <-
        self$options$supportSignalVectorN > 1L
      
      reestimationPlotsAvailable <-
        self$options$twostepsN > 0L

      parseNumericVector <- function(x) {
        parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])

        if (length(parts) == 0 || any(parts == "")) {
          return(numeric(0))
        }

        suppressWarnings(as.numeric(parts))
      }

      Mvalues <- unique(parseNumericVector(self$options$M))
      cvPlotAvailable <- length(Mvalues) > 1L
      
      self$results$plot1$setVisible(isTRUE(self$options$plot1))
      self$results$plot2$setVisible(
        isTRUE(self$options$plot2) && supportGridAvailable
      )
      self$results$plot3$setVisible(
        isTRUE(self$options$plot3) && supportGridAvailable
      )
      self$results$plot4$setVisible(
        isTRUE(self$options$plot4) && supportGridAvailable
      )
      self$results$plot5$setVisible(
        isTRUE(self$options$plot5) && supportGridAvailable
      )
      self$results$plot6$setVisible(
        isTRUE(self$options$plot6) && reestimationPlotsAvailable
      )
      self$results$plot7$setVisible(
        isTRUE(self$options$plot7) && reestimationPlotsAvailable
      )

      needsTrueCoef <- (isTRUE(self$options$plot5) && supportGridAvailable) ||
        (isTRUE(self$options$plot7) && reestimationPlotsAvailable)
      
      trueCoef <- NULL
      
      vars <- c(dep, covs, factors)
      
      data <- self$data[, vars, drop = FALSE]
      data <- data[stats::complete.cases(data), , drop = FALSE]
      
      if (nrow(data) < 3) {
        jmvcore::reject("This analysis requires at least 3 complete cases ({n} found).",
          code = "too_few_cases", n = nrow(data))
      }
      
      for (f in factors) {
        data[[f]] <- droplevels(as.factor(data[[f]]))
      }
      
      constantCovs <- covs[
        vapply(
          data[, covs, drop = FALSE],
          function(x) length(unique(x)) < 2,
          logical(1)
        )
      ]
      
      constantFactors <- factors[
        vapply(
          data[, factors, drop = FALSE],
          function(x) nlevels(x) < 2,
          logical(1)
        )
      ]
      
      constantVars <- c(constantCovs, constantFactors)
      
      if (length(constantVars) > 0) {
        jmvcore::reject("The following predictors have no variation: {variables}.",
          code = "constant_predictors", variables = paste(constantVars, collapse = ", "))
      }
      
      terms <- c(covs, factors)
      
      formStr <- jmvcore::constructFormula(
        dep = dep,
        terms = c(covs, factors)
      )
      
      if (!isTRUE(self$options$intercept))
        formStr <- paste(formStr, "- 1")
      
      form <- stats::as.formula(formStr)
      
      precisionWarning <- NULL
      if (needsTrueCoef) {
        validated <- .gceTrueCoefficients(self$options$trueCoef, form, data)
        trueCoef <- validated$values
        precisionWarning <- validated$warning
      }

      if (self$options$supportSignalVectorMin >= self$options$supportSignalVectorMax) {
        jmvcore::reject("The minimum flexibility parameter must be smaller than the maximum.", code = "invalid_signal_range")
      }
      
      if (self$options$cvNfolds > nrow(data)) {
        jmvcore::reject("The number of CV folds cannot exceed the number of complete cases ({n}).",
          code = "too_many_folds", n = nrow(data))
      }
      
      Jvalues <- unique(parseNumericVector(self$options$J))
      weightValues <- unique(parseNumericVector(self$options$weight))
      
      self$results$plotCV$setVisible(
        isTRUE(self$options$plotCV) && cvPlotAvailable
      )
      
      if (
        length(Mvalues) == 0 ||
        anyNA(Mvalues) ||
        any(!is.finite(Mvalues)) ||
        any(Mvalues < 3) ||
        any(Mvalues != floor(Mvalues)) ||
        any(Mvalues %% 2 == 0)
      ) {
        jmvcore::reject("Signal support points must be comma-separated odd integers greater than or equal to 3.", code = "invalid_signal_points")
      }
      
      if (
        length(Jvalues) == 0 ||
        anyNA(Jvalues) ||
        any(!is.finite(Jvalues)) ||
        any(Jvalues < 3) ||
        any(Jvalues != floor(Jvalues)) ||
        any(Jvalues %% 2 == 0)
      ) {
        jmvcore::reject("Noise support points must be comma-separated odd integers greater than or equal to 3.", code = "invalid_noise_points")
      }
      
      if (
        length(weightValues) == 0 ||
        anyNA(weightValues) ||
        any(!is.finite(weightValues)) ||
        any(weightValues < 0) ||
        any(weightValues > 1)
      ) {
        jmvcore::reject("Noise weight values must be comma-separated numbers between 0 and 1.", code = "invalid_noise_weights")
      }
      
      bootB <- 0
      bootMethod <- "residuals"
      
      if (self$options$bootstrap) {
        
        bootB <- self$options$bootB
        bootMethod <- self$options$bootMethod
        
      }
      
      fitCache <- self$results$fitCache$state
      cached <- fitCache
      
      useCachedFit <- is.list(cached) &&
        identical(cached$version, 4L) &&
        is.list(cached$fit$best) &&
        is.matrix(cached$coefficients) &&
        !is.null(cached$ci)

      self$results$text$setContent(if (useCachedFit) {
        "<p>Updating results from the fitted model...</p>"
      } else {
        "<p>Starting fit...</p>"
      })
      private$.checkpoint()

      start_time <- Sys.time()
      
      warningMsg <- NULL
      
      fixedSupportSignal <- NULL
      
      if (self$options$supportSignalVectorN == 1L) {
        
        fixedSupportSignal <- self$options$supportSignal
        
        if (
          length(fixedSupportSignal) != 1L ||
          is.na(fixedSupportSignal) ||
          !is.finite(fixedSupportSignal) ||
          fixedSupportSignal <= 0
        ) {
          jmvcore::reject("The flexibility parameter must be a positive finite number.",
            code = "invalid_signal_support")
        }
      }
      
      if (useCachedFit) {
        fit <- cached$fit
        warningMsg <- fitCache$warningMsg
        fitElapsed <- fitCache$fitElapsed
      } else {
        fit <- tryCatch(
          withCallingHandlers(
          GCEstim::cv.lmgce(
          formula = form,
          data = data,
          support.method = "standardized",
          support.signal = fixedSupportSignal,
          support.signal.vector.n = self$options$supportSignalVectorN,
          support.signal.vector.min = self$options$supportSignalVectorMin,
          support.signal.vector.max = self$options$supportSignalVectorMax,
          support.signal.points = Mvalues,
          support.noise.points = Jvalues,
          cv.nfolds = self$options$cvNfolds,
          errormeasure = self$options$errorMeasure,
          errormeasure.which = self$options$errorMeasureWhich,
          weight = weightValues,
          twosteps.n = self$options$twostepsN,
          method = self$options$method,
          caseGLM = "D",
          boot.B = bootB,
          boot.method = bootMethod,
          OLS = TRUE,
          seed = self$options$seed
        ),
        warning = function(w) {
          message <- conditionMessage(w)
          
          warningMsg <<- if (is.null(warningMsg) || !nzchar(warningMsg)) {
            message
          } else {
            paste(warningMsg, message, sep = "; ")
          }
          
          invokeRestart("muffleWarning")
        }
        ),
          error = function(e) e
        )

        if (inherits(fit, "error")) {
          jmvcore::reject("Model estimation failed: {message}",
            code = "fit_failed", message = conditionMessage(fit))
        }
      }

      fitWarningMsg <- warningMsg

      bestFit <- fit$best

      if (is.null(bestFit)) {
        jmvcore::reject("Cross-validation did not return a selected model.", code = "no_selected_model")
      }
      
      compactFit <- if (useCachedFit) fit else .gceCompactCVFit(fit)

      if (!useCachedFit)
        fitElapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

      supportLimits <- as.matrix(bestFit$support.matrix)
      supportTable <- self$results$signalSupports
      supportTable$setTitle(sprintf(
        "Signal Supports (%d equally spaced points)",
        as.integer(fit$support.signal.points.best)
      ))
      if ((length(supportTable$rowKeys) == 0L || !supportTable$isFilled())) for (i in seq_len(nrow(supportLimits))) {
        supportTable$addRow(rowKey = as.character(i), values = list(
          term = rownames(supportLimits)[i],
          lower = as.numeric(supportLimits[i, 1]),
          upper = as.numeric(supportLimits[i, 2])
        ))
      }

      if (self$options$saveFitted &&
          self$results$saveFitted$isNotFilled()) {
        self$results$saveFitted$setRowNums(rownames(data))
        self$results$saveFitted$setValues(
          as.numeric(bestFit$fitted.values)
        )
      }
      
      if (self$options$saveResiduals &&
          self$results$saveResiduals$isNotFilled()) {
        self$results$saveResiduals$setRowNums(rownames(data))
        self$results$saveResiduals$setValues(
          as.numeric(bestFit$residuals)
        )
      }
      
      if (self$options$saveW) {
        w <- as.matrix(bestFit$w)
        
        if (nrow(w) != nrow(data) || ncol(w) < 1L) {
          message <- "Noise probabilities could not be saved because their dimensions do not match the complete cases."
          warningMsg <- if (is.null(warningMsg) || !nzchar(warningMsg)) {
            message
          } else {
            paste(warningMsg, message, sep = "; ")
          }
        } else if (self$results$saveW$isNotFilled()) {
          keys <- as.character(seq_len(ncol(w)))
          
          self$results$saveW$set(
            keys,
            paste0("TASW-GCE w", keys),
            paste0("Selected model noise probability at support point ", keys),
            rep("continuous", length(keys))
          )
          
          self$results$saveW$setRowNums(rownames(data))
          
          for (j in seq_len(ncol(w))) {
            self$results$saveW$setValues(
              as.numeric(w[, j]),
              index = j
            )
          }
        }
      }
      
      elapsed <- fitElapsed
      
      time_txt <- if (elapsed < 60)
        paste0(round(elapsed, 3), " s")
      else
        paste0(
          floor(elapsed / 60), " min ",
          round(elapsed %% 60, 1), " s"
        )
      
      self$results$text$setContent("<p>Fit completed. Computing CI...</p>")
      private$.checkpoint()
      
      ## Confidence Interval
      
      if (useCachedFit) {
        ci <- cached$ci
      } else {
        if (self$options$bootstrap) {
          ci <- tryCatch(
            stats::confint(
              bestFit,
              level = self$options$bootConfLevel,
              method = self$options$bootCIMethod
            ),
            error = function(e) e
          )
        } else {
          ci <- tryCatch(
            stats::confint(
              bestFit,
              level = self$options$bootConfLevel
            ),
            error = function(e) e
          )
        }

        if (inherits(ci, "error")) {
          jmvcore::reject("Confidence intervals: {message}",
            code = "ci_failed", message = conditionMessage(ci))
        }


      }
      
      self$results$text$setContent("<p>CI completed. Filling tables...</p>")
      private$.checkpoint()
      
      ## Summary table
      
      summ <- if (useCachedFit) {
        list(coefficients = cached$coefficients)
      } else {
        tryCatch(
          get("summary.lmgce", envir = asNamespace("GCEstim"),
              inherits = FALSE)(bestFit),
          error = function(e) e
        )
      }

      if (inherits(summ, "error")) {
        jmvcore::reject("Model summary: {message}",
          code = "summary_failed", message = conditionMessage(summ))
      }
      
      summaryTable <- self$results$modelSummary
      
      summaryTable$setRow(rowKey = "formula", values = list(
        measure = "Formula",
        value = paste(deparse(form), collapse = " ")
      ))
      
      summaryTable$setRow(rowKey = "n", values = list(
        measure = "Complete cases",
        value = as.character(nrow(data))
      ))
      
      summaryTable$setRow(rowKey = "cvNfolds", values = list(
        measure = "Cross-validation folds",
        value = as.character(self$options$cvNfolds)
      ))
      
      summaryTable$setRow(rowKey = "errorMeasure", values = list(
        measure = "Prediction-error measure",
        value = self$options$errorMeasure
      ))
      
      summaryTable$setRow(rowKey = "errorMeasureWhich", values = list(
        measure = "Selection rule",
        value = c(
          "1se" = "Within one standard error",
          "min" = "Minimum error"
        )[[self$options$errorMeasureWhich]]
      ))
      
      summaryTable$setRow(rowKey = "seed", values = list(
        measure = "Seed",
        value = as.character(self$options$seed)
      ))
      
      methodLabel <- sub("^dual\\.", "Dual: ", self$options$method)
      methodLabel <- sub("^primal\\.", "Primal: ", methodLabel)
      
      summaryTable$setRow(rowKey = "method", values = list(
        measure = "Optimization method",
        value = methodLabel
      ))
      
      summaryTable$setRow(rowKey = "supportMethod", values = list(
        measure = "Support method",
        value = "Standardized coefficients"
      ))
      
      summaryTable$setRow(rowKey = "signalSupportSpecification", values = list(
        measure = "Signal support specification",
        value = "Symmetric"
      ))
      
      noiseResponse <- as.numeric(data[[dep]])
      
      noiseLimit <- 3 * stats::sd(noiseResponse)
      
      noiseInterval <- paste0(
        "[",
        format(signif(-noiseLimit, 6), trim = TRUE), ", ",
        format(signif(noiseLimit, 6), trim = TRUE), "]"
      )
      
      summaryTable$setRow(rowKey = "noiseSupportMethod", values = list(
        measure = "Noise support specification",
        value = paste0("3 sigma: ", noiseInterval)
      ))
      
      summaryTable$setRow(rowKey = "supportSignalVectorN", values = list(
        measure = "Number of support spaces",
        value = as.character(self$options$supportSignalVectorN)
      ))
      
      summaryTable$setRow(rowKey = "M", values = list(
        measure = "Candidate signal support points",
        value = paste(Mvalues, collapse = ", ")
      ))
      
      summaryTable$setRow(rowKey = "J", values = list(
        measure = "Candidate noise support points",
        value = paste(Jvalues, collapse = ", ")
      ))
      
      summaryTable$setRow(rowKey = "weight", values = list(
        measure = "Candidate noise weight values",
        value = paste(weightValues, collapse = ", ")
      ))
      
      summaryTable$setRow(rowKey = "twostepsN", values = list(
        measure = "Post-GCE re-estimations",
        value = as.character(self$options$twostepsN)
      ))
      
      summaryTable$setRow(rowKey = "ciType", values = list(
        measure = "Confidence interval type",
        value = if (self$options$bootstrap) {
          c(
            percentile = "Bootstrap percentile",
            basic = "Bootstrap basic"
          )[[self$options$bootCIMethod]]
        } else {
          "Asymptotic normal"
        }
      ))
      
      summaryTable$setRow(rowKey = "ciLevel", values = list(
        measure = "Confidence level",
        value = paste0(
          100 * self$options$bootConfLevel,
          "%"
        )
      ))
      
      if (self$options$bootstrap) {
        
        summaryTable$setRow(rowKey = "bootB", values = list(
          measure = "Bootstrap samples",
          value = as.character(self$options$bootB)
        ))
        
        summaryTable$setRow(rowKey = "bootMethod", values = list(
          measure = "Bootstrap resampling method",
          value = c(
            residuals = "Residual",
            cases = "Case",
            wild = "Wild"
          )[[self$options$bootMethod]]
        ))
      }
      
      ## Cross-validation results table
      
      cvResults <- as.data.frame(fit$results)
      
      cvResults <- cvResults[
        order(cvResults$error.measure.cv.mean, na.last = TRUE),
        ,
        drop = FALSE
      ]
      
      nCVRows <- min(
        self$options$cvTableRows,
        nrow(cvResults)
      )
      
      measure <- as.character(self$options$errorMeasure)
      cvTable <- self$results$cvResults
      
      cvTable$getColumn("errorMeasure")$setTitle(paste("In-sample", measure))
      cvTable$getColumn("cvMean")$setTitle(paste("CV", measure))
      
      if (nCVRows > 0) {
        
        for (i in seq_len(nCVRows)) {

          if (cvTable$isFilled(col = "cvMean", rowKey = as.character(i)))
            next
          
          cvTable$setRow(
            rowKey = as.character(i),
            values = list(
              signalPoints = cvResults$support.signal.points[i],
              noisePoints = cvResults$support.noise.points[i],
              alpha = cvResults$weight[i],
              errorMeasure = cvResults$error.measure[i],
              cvMean = cvResults$error.measure.cv.mean[i],
              convergence = if (cvResults$convergence[i] == 0) "Yes" else "No",
              time = cvResults$time[i]
            )
          )
        }
      }
      
      ## Coefficients table

      coefs <- summ$coefficients
      coefTable <- self$results$coefficients
      
      if ((length(coefTable$rowKeys) == 0L || !coefTable$isFilled())) for (i in seq_len(nrow(coefs))) {
        
        term <- rownames(coefs)[i]
        
        coefTable$addRow(
          rowKey = term,
          values = list(
            term = term,
            estimate = coefs[i, "Estimate"],
            ciLower = ci[term, 1],
            ciUpper = ci[term, 2],
            se = coefs[i, "Std. Deviation"],
            z = coefs[i, "z value"],
            p = coefs[i, "Pr(>|t|)"]
          )
        )
      }
      
      for (key in grep("^p[0-9]+$", names(self$results$signalProbabilities$columns), value = TRUE))
        self$results$signalProbabilities$getColumn(key)$setTitle(
          paste0("p (point ", substring(key, 2L), ")"))
      pMatrix <- as.matrix(bestFit$p)
      
      if (!is.null(bestFit$p) &&
          nrow(pMatrix) == nrow(coefs) &&
          ncol(pMatrix) > 0L &&
          (length(self$results$signalProbabilities$rowKeys) == 0L ||
           !self$results$signalProbabilities$isFilled())) {
        
        pTable <- self$results$signalProbabilities
        pKeys <- paste0("p", seq_len(ncol(pMatrix)))
        
        if (length(pTable$rowKeys) == 0L) for (j in seq_len(ncol(pMatrix))) {
          pTable$addColumn(
            name = pKeys[j],
            title = paste0("p (point ", j, ")"),
            type = "number",
            format = "zto"
          )
        }
        
        pTerms <- rownames(pMatrix)
        if (is.null(pTerms))
          pTerms <- rownames(coefs)
        
        for (i in seq_len(nrow(pMatrix))) {
          pTable$addRow(
            rowKey = as.character(i),
            values = c(
              list(term = pTerms[i]),
              stats::setNames(
                as.list(as.numeric(pMatrix[i, ])),
                pKeys
              )
            )
          )
        }
      }
      
      .gceLabelTables(self$results, rownames(coefs),
        supportTerms = rownames(supportLimits),
        probabilityTerms = rownames(pMatrix) %||% rownames(coefs))

      self$results$text$setContent("<p>Tables completed. Plotting...</p>")
      private$.checkpoint()
      
      if (useCachedFit) {
        plotStates <- cached$plots
        precisionBasis <- cached$precisionBasis
        plotWarnings <- fitCache$plotWarnings %||% character()
      } else {
        prepared <- .gceCachePlotStates(
          fit, ci, self$options$bootConfLevel,
          if (self$options$bootstrap) self$options$bootCIMethod else "z",
          supportGridAvailable, reestimationPlotsAvailable,
          cvPlotAvailable
        )
        plotStates <- prepared$states
        precisionBasis <- .gcePrecisionBasis(bestFit)
        plotWarnings <- prepared$warnings
      }

      # Canonical plot data lives once in fitCache. Images retain references.
      .gceEntropyWarning(self$results, self$options, plotStates)
      for (key in c("plot1", "plot2", "plot3", "plot4", "plot6", "plotCV")) {
        if (isTRUE(self$options[[key]]) && !is.null(plotStates[[key]]))
          self$results[[key]]$setState(list(key = key,
            OLS = isTRUE(self$options$OLS), NormEnt = isTRUE(self$options$NormEnt)))
      }
      for (which in c(5L, 7L)) {
        key <- paste0("plot", which)
        available <- if (which == 5L) supportGridAvailable else reestimationPlotsAvailable
        if (!isTRUE(self$options[[key]]) || !available) next
        state <- NULL
        if (!is.null(trueCoef)) {
          state <- tryCatch(.gcePrecisionState(precisionBasis, which, trueCoef),
            error = function(e) {
              precisionWarning <<- paste0("Precision-error plot skipped: ", conditionMessage(e))
              NULL
            })
          if (!is.null(state)) {
            if (!isTRUE(self$options$OLS))
              state$overlays <- Filter(function(x) !"yintercept" %in% names(x), state$overlays)
            state$showNEP <- isTRUE(self$options$NormEnt)
          }
        }
        self$results[[key]]$setState(state)
      }

      if (!useCachedFit) {
        self$results$fitCache$setState(list(
          version = 4L,
          fit = compactFit, ci = ci,
          coefficients = as.matrix(summ$coefficients),
          plots = plotStates, precisionBasis = precisionBasis,
          warningMsg = fitWarningMsg,
          plotWarnings = plotWarnings,
          fitElapsed = fitElapsed
        ))
      }

      ## Complete
      
      warningMessages <- c(plotWarnings, precisionWarning)
      
      if (!is.null(warningMsg) && nzchar(warningMsg))
        warningMessages <- c(warningMsg, warningMessages)
      
      if (length(warningMessages)) {
        warningNotice <- jmvcore::Notice$new(
          options = self$options,
          name = ".fitWarning",
          type = jmvcore::NoticeType$WARNING
        )
        warningNotice$setContent(
          paste(unique(warningMessages), collapse = "; ")
        )
        self$results$insert(1, warningNotice)
      }
      
      timeNotice <- jmvcore::Notice$new(
        options = self$options,
        name = ".fitTime",
        type = jmvcore::NoticeType$INFO
      )
      timeNotice$setContent(if (useCachedFit) {
        paste0("Previous fit reused (original fit time: ", time_txt, ").")
      } else {
        paste0("Computation time: ", time_txt)
      })
      self$results$insert(1, timeNotice)
      
      self$results$text$setVisible(FALSE)
      
    },
    
    .plot1 = function(image, ggtheme, theme)
      .gceRenderEstimates(.gcePlotImage(image, self$results$fitCache$state), ggtheme, theme),
    .plot2 = function(image, ggtheme, theme)
      .gceRenderDiagnostic(.gcePlotImage(image, self$results$fitCache$state), ggtheme, theme),
    .plot3 = function(image, ggtheme, theme)
      .gceRenderDiagnostic(.gcePlotImage(image, self$results$fitCache$state), ggtheme, theme),
    .plot4 = function(image, ggtheme, theme)
      .gceRenderDiagnostic(.gcePlotImage(image, self$results$fitCache$state), ggtheme, theme),
    .plot5 = function(image, ggtheme, theme)
      .gceRenderDiagnostic(.gcePlotImage(image, self$results$fitCache$state), ggtheme, theme),
    .plot6 = function(image, ggtheme, theme)
      .gceRenderDiagnostic(.gcePlotImage(image, self$results$fitCache$state), ggtheme, theme),
    .plot7 = function(image, ggtheme, theme)
      .gceRenderDiagnostic(.gcePlotImage(image, self$results$fitCache$state), ggtheme, theme),
    .plotCV = function(image, ggtheme, theme)
      .gceRenderCV(.gcePlotImage(image, self$results$fitCache$state), ggtheme, theme)
  )
)
