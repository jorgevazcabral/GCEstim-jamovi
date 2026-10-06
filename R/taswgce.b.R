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
          note = "Asymptotic Std. Deviation, z, and p are hidden when bootstrap intervals are selected.",
          init = TRUE
        )
      }
    },
    
    .run = function() {
        
      dep  <- unlist(self$options$dep)
      covs    <- unlist(self$options$covs)
      factors <- unlist(self$options$factors)
      
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
      
      selectedMainPlots <- c(
        self$options$plot1,
        self$options$plot2 && supportGridAvailable,
        self$options$plot3 && supportGridAvailable,
        self$options$plot4 && supportGridAvailable,
        self$options$plot5 && supportGridAvailable,
        self$options$plot6 && reestimationPlotsAvailable,
        self$options$plot7 && reestimationPlotsAvailable
      )

      combineMainPlots <-
        sum(selectedMainPlots) +
        as.integer(isTRUE(self$options$plotCV) && cvPlotAvailable) > 1L

      self$results$plotCombined$setVisible(combineMainPlots)
      self$results$plot1$setVisible(
        isTRUE(self$options$plot1) && !combineMainPlots
      )
      
      self$results$plot2$setVisible(
        isTRUE(self$options$plot2) && supportGridAvailable && !combineMainPlots
      )
      
      self$results$plot3$setVisible(
        isTRUE(self$options$plot3) && supportGridAvailable && !combineMainPlots
      )
      
      self$results$plot4$setVisible(
        isTRUE(self$options$plot4) && supportGridAvailable && !combineMainPlots
      )
      
      self$results$plot5$setVisible(
        isTRUE(self$options$plot5) && supportGridAvailable && !combineMainPlots
      )
      
      self$results$plot6$setVisible(
        isTRUE(self$options$plot6) && reestimationPlotsAvailable && !combineMainPlots
      )
      
      self$results$plot7$setVisible(
        isTRUE(self$options$plot7) && reestimationPlotsAvailable && !combineMainPlots
      )
      
      needsTrueCoef <- self$options$plot5 || self$options$plot7
      
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
      
      form <- stats::reformulate(
        termlabels = terms,
        response = dep,
        intercept = self$options$intercept
      )
      
      if (needsTrueCoef) {
        
        trueCoefText <- trimws(self$options$trueCoef)
        
        if (trueCoefText == "") {
          jmvcore::reject("Plots 5 and 7 require a vector of true coefficients.", code = "true_coef_required")
        }
        
        trueCoef <- suppressWarnings(
          as.numeric(strsplit(trueCoefText, ",")[[1]])
        )
        
        if (any(is.na(trueCoef))) {
          jmvcore::reject("True coefficients must be numeric and separated by commas.", code = "invalid_true_coef")
        }
        
        X <- stats::model.matrix(form, data)
        
        expectedCoef <- ncol(X)
        
        if (length(trueCoef) != expectedCoef) {
          jmvcore::reject("The true coefficient vector must contain {n} values (including the intercept).",
            code = "true_coef_length", n = expectedCoef)
        }
      }
      
      if (self$options$supportSignalVectorMin >= self$options$supportSignalVectorMax) {
        jmvcore::reject("The minimum signal support range must be smaller than the maximum.", code = "invalid_signal_range")
      }
      
      if (self$options$cvNfolds > nrow(data)) {
        jmvcore::reject("The number of CV folds cannot exceed the number of complete cases ({n}).",
          code = "too_many_folds", n = nrow(data))
      }
      
      Jvalues <- unique(parseNumericVector(self$options$J))
      weightValues <- unique(parseNumericVector(self$options$weight))
      
      self$results$plotCV$setVisible(
        isTRUE(self$options$plotCV) && cvPlotAvailable && !combineMainPlots
      )
      
      self$results$plotCV$setSize(
        400,
        160 + 140 * length(weightValues)
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
      
      self$results$text$setContent("<p>Starting fit...</p>")
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
          jmvcore::reject("The signal support value must be a positive finite number.",
            code = "invalid_signal_support")
        }
      }
      
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
          caseGLM = self$options$caseGLM,
          boot.B = bootB,
          boot.method = bootMethod,
          OLS = self$options$OLS,
          seed = self$options$seed
        ),
        warning = function(w) {
          warningMsg <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        }
        ),
        error = function(e) e
      )
      
      if (inherits(fit, "error")) {
        jmvcore::reject("Model estimation failed: {message}",
          code = "fit_failed", message = conditionMessage(fit))
      }
      
      bestFit <- fit$best
      
      if (is.null(bestFit)) {
        jmvcore::reject("Cross-validation did not return a selected model.", code = "no_selected_model")
      }
      
      supportLimits <- as.matrix(bestFit$support.matrix)
      supportTable <- self$results$signalSupports
      supportTable$setTitle(sprintf(
        "Signal Supports (%d equally spaced points)",
        as.integer(fit$support.signal.points.best)
      ))
      for (i in seq_len(nrow(supportLimits))) {
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
      
      if (self$options$saveW && self$options$caseGLM == "D") {
        w <- as.matrix(bestFit$w)
        
        if (nrow(w) != nrow(data) || ncol(w) < 1L) {
          jmvcore::reject("Noise probabilities do not match the complete cases.", code = "noise_probability_dimensions")
        }
        
        if (self$results$saveW$isNotFilled()) {
          keys <- as.character(seq_len(ncol(w)))
          
          self$results$saveW$set(
            keys,
            paste0("TARW-GCE w", keys),
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
      
      elapsed <- as.numeric(
        difftime(Sys.time(), start_time, units = "secs")
      )
      
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
      
      ci <- NULL
      
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
      
      self$results$text$setContent("<p>CI completed. Filling tables...</p>")
      private$.checkpoint()
      
      ## Summary table
      
      summ <- tryCatch(
        summary(bestFit),
        error = function(e) e
      )
      
      if (inherits(summ, "error")) {
        jmvcore::reject("Model summary: {message}",
          code = "summary_failed", message = conditionMessage(summ))
      }
      
      summaryTable <- self$results$modelSummary
      
      summaryTable$addRow(rowKey = "formula", values = list(
        measure = "Formula",
        value = paste(deparse(form), collapse = " ")
      ))
      
      summaryTable$addRow(rowKey = "n", values = list(
        measure = "Complete cases",
        value = as.character(nrow(data))
      ))
      
      summaryTable$addRow(rowKey = "cvNfolds", values = list(
        measure = "Cross-validation folds",
        value = as.character(self$options$cvNfolds)
      ))
      
      summaryTable$addRow(rowKey = "errorMeasure", values = list(
        measure = "Prediction-error measure",
        value = self$options$errorMeasure
      ))
      
      summaryTable$addRow(rowKey = "errorMeasureWhich", values = list(
        measure = "Selection rule",
        value = self$options$errorMeasureWhich
      ))
      
      summaryTable$addRow(rowKey = "seed", values = list(
        measure = "Seed",
        value = as.character(self$options$seed)
      ))
      
      summaryTable$addRow(rowKey = "method", values = list(
        measure = "Method",
        value = self$options$method
      ))
      
      summaryTable$addRow(rowKey = "supportMethod", values = list(
        measure = "Support method",
        value = "Standardized coefficients"
      ))
      
      summaryTable$addRow(rowKey = "signalSupportSpecification", values = list(
        measure = "Signal support specification",
        value = "Symmetric"
      ))
      
      noiseResponse <- as.numeric(data[[dep]])
      if (self$options$caseGLM %in% c("M", "NM")) {
        noiseX <- stats::model.matrix(form, data)
        noiseResponse <- as.numeric(crossprod(noiseX, noiseResponse))
        if (self$options$caseGLM == "NM")
          noiseResponse <- noiseResponse / nrow(noiseX)
      }
      
      noiseLimit <- 3 * stats::sd(noiseResponse)
      
      noiseInterval <- paste0(
        "[",
        format(signif(-noiseLimit, 6), trim = TRUE), ", ",
        format(signif(noiseLimit, 6), trim = TRUE), "]"
      )
      
      summaryTable$addRow(rowKey = "noiseSupportMethod", values = list(
        measure = "Noise support specification",
        value = paste0("3 sigma: ", noiseInterval)
      ))
      
      summaryTable$addRow(rowKey = "supportSignalVectorN", values = list(
        measure = "Number of support spaces",
        value = as.character(self$options$supportSignalVectorN)
      ))
      
      summaryTable$addRow(rowKey = "M", values = list(
        measure = "Candidate signal support points",
        value = paste(Mvalues, collapse = ", ")
      ))
      
      summaryTable$addRow(rowKey = "J", values = list(
        measure = "Candidate noise support points",
        value = paste(Jvalues, collapse = ", ")
      ))
      
      summaryTable$addRow(rowKey = "weight", values = list(
        measure = "Candidate noise weight values",
        value = paste(weightValues, collapse = ", ")
      ))
      
      summaryTable$addRow(rowKey = "caseGLM", values = list(
        measure = "GLM case",
        value = self$options$caseGLM
      ))
      
      summaryTable$addRow(rowKey = "twostepsN", values = list(
        measure = "Post-GCE re-estimations",
        value = as.character(self$options$twostepsN)
      ))
      
      summaryTable$addRow(rowKey = "ciType", values = list(
        measure = "Confidence interval type",
        value = if (self$options$bootstrap) {
          "Bootstrap"
        } else {
          "Asymptotic normal"
        }
      ))
      
      summaryTable$addRow(rowKey = "ciMethod", values = list(
        measure = "Confidence interval method",
        value = if (self$options$bootstrap) {
          self$options$bootCIMethod
        } else {
          "z"
        }
      ))
      
      summaryTable$addRow(rowKey = "ciLevel", values = list(
        measure = "Confidence level",
        value = paste0(
          100 * self$options$bootConfLevel,
          "%"
        )
      ))
      
      if (self$options$bootstrap) {
        
        summaryTable$addRow(rowKey = "bootB", values = list(
          measure = "Bootstrap samples",
          value = as.character(self$options$bootB)
        ))
        
        summaryTable$addRow(rowKey = "bootMethod", values = list(
          measure = "Bootstrap resampling method",
          value = self$options$bootMethod
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
      
      cvTable <- self$results$cvResults
      
      if (nCVRows > 0) {
        
        for (i in seq_len(nCVRows)) {
          
          cvTable$addRow(
            rowKey = paste0("cv_", i),
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
      
      hasCI <- self$options$bootstrap && !is.null(ci)
      
      for (i in seq_len(nrow(coefs))) {
        
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
      
      pMatrix <- as.matrix(bestFit$p)
      
      if (!is.null(bestFit$p) &&
          nrow(pMatrix) == nrow(coefs) &&
          ncol(pMatrix) > 0L) {
        
        pTable <- self$results$signalProbabilities
        pKeys <- paste0("p", seq_len(ncol(pMatrix)))
        
        for (j in seq_len(ncol(pMatrix))) {
          pTable$addColumn(
            name = pKeys[j],
            title = paste0("p_", j),
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
      
      self$results$text$setContent("<p>Tables completed. Plotting...</p>")
      private$.checkpoint()
      
      plots <- list()

      if (self$options$plot1)
        plots$plot1 <- tryCatch(plot(bestFit,
                                     which = 1,
                                     ci.level = self$options$bootConfLevel,
                                     ci.method = ifelse(self$options$bootstrap,
                                                        self$options$bootCIMethod,
                                                        "z"),
                                     OLS = self$options$OLS),
                                error = function(e) NULL)

      if (self$options$plot2 && supportGridAvailable)
        plots$plot2 <- tryCatch(plot(bestFit,
                                     which = 2,
                                     OLS = self$options$OLS,
                                     NormEnt = self$options$NormEnt),
                                error = function(e) NULL)

      if (self$options$plot3 && supportGridAvailable)
        plots$plot3 <- tryCatch(plot(bestFit,
                                     which = 3,
                                     OLS = self$options$OLS),
                                error = function(e) NULL)

      if (self$options$plot4 && supportGridAvailable)
        plots$plot4 <- tryCatch(plot(bestFit,
                                     which = 4,
                                     OLS = self$options$OLS),
                                error = function(e) NULL)

      if (self$options$plot5 && supportGridAvailable)
        plots$plot5 <- tryCatch(plot(bestFit,
                                     which = 5,
                                     coef = trueCoef,
                                     OLS = self$options$OLS,
                                     NormEnt = self$options$NormEnt),
                                error = function(e) NULL)

      if (self$options$plot6 && reestimationPlotsAvailable)
        plots$plot6 <- tryCatch(plot(bestFit,
                                     which = 6,
                                     OLS = self$options$OLS),
                                error = function(e) NULL)

      if (self$options$plot7 && reestimationPlotsAvailable)
        plots$plot7 <- tryCatch(plot(bestFit,
                                     which = 7,
                                     coef = trueCoef,
                                     OLS = self$options$OLS),
                                error = function(e) NULL)
      
      if (isTRUE(self$options$plotCV) && cvPlotAvailable) {
        plots$plotCV <- tryCatch(
          plot(
            fit,
            which = 1,
            ncol = 1,
            scales = "free"
          ),
          error = function(e) NULL
        )
      }

      if (combineMainPlots) {
        plotNames <- c(
          if (isTRUE(self$options$plotCV) && cvPlotAvailable) "plotCV",
          paste0("plot", which(selectedMainPlots))
        )
        plotItems <- lapply(plots[plotNames], function(x) {
          if (inherits(x, "ggplot"))
            return(x)
          if (is.null(x) || length(x) == 0L ||
              !inherits(x[[1L]], "ggplot"))
            return(NULL)
          x[[1L]]
        })
        plotItems <- Filter(Negate(is.null), plotItems)

        if (length(plotItems) > 0L) {
          plotHeights <- ifelse(
            names(plotItems) == "plotCV",
            160 + 140 * length(weightValues),
            320
          )
          self$results$plotCombined$setSize(600, sum(plotHeights))
          self$results$plotCombined$setState(list(
            p = ggpubr::ggarrange(
              plotlist = unname(plotItems), ncol = 1,
              heights = unname(plotHeights)
            )
          ))
        } else {
          self$results$plotCombined$setState(list(p = NULL))
        }
      }
      
      if (self$options$plot1 && !combineMainPlots)
        self$results$plot1$setState(list(p = plots$plot1))

      if (self$options$plot2 && supportGridAvailable && !combineMainPlots)
        self$results$plot2$setState(list(p = plots$plot2))
      
      if (self$options$plot3 && supportGridAvailable && !combineMainPlots)
        self$results$plot3$setState(list(p = plots$plot3))
      
      if (self$options$plot4 && supportGridAvailable && !combineMainPlots)
        self$results$plot4$setState(list(p = plots$plot4))
      
      if (self$options$plot5 && supportGridAvailable && !combineMainPlots)
        self$results$plot5$setState(list(p = plots$plot5))

      if (self$options$plot6 && reestimationPlotsAvailable && !combineMainPlots)
        self$results$plot6$setState(list(p = plots$plot6))

      if (self$options$plot7 && reestimationPlotsAvailable && !combineMainPlots)
        self$results$plot7$setState(list(p = plots$plot7))
      
      if (isTRUE(self$options$plotCV) && cvPlotAvailable && !combineMainPlots)
        self$results$plotCV$setState(list(p = plots$plotCV))
      
      ## Complete
      
      if (!is.null(warningMsg) && nzchar(warningMsg)) {
        warningNotice <- jmvcore::Notice$new(
          options = self$options,
          name = ".fitWarning",
          type = jmvcore::NoticeType$WARNING
        )
        warningNotice$setContent(warningMsg)
        self$results$insert(1, warningNotice)
      }
      
      timeNotice <- jmvcore::Notice$new(
        options = self$options,
        name = ".fitTime",
        type = jmvcore::NoticeType$INFO
      )
      timeNotice$setContent(paste0("Computation time: ", time_txt))
      self$results$insert(1, timeNotice)
      
      self$results$text$setVisible(FALSE)
      
    },
    
    .plot1 = function(image, ggtheme, theme) private$.printPlot(image),
    .plot2 = function(image, ggtheme, theme) private$.printPlot(image),
    .plot3 = function(image, ggtheme, theme) private$.printPlot(image),
    .plot4 = function(image, ggtheme, theme) private$.printPlot(image),
    .plot5 = function(image, ggtheme, theme) private$.printPlot(image),
    .plot6 = function(image, ggtheme, theme) private$.printPlot(image),
    .plot7 = function(image, ggtheme, theme) private$.printPlot(image),
    .plotCombined = function(image, ggtheme, theme) private$.printPlot(image),
    .plotCV = function(image, ggtheme, theme) private$.printPlot(image),

    .printPlot = function(image) {

      p <- image$state$p

      if (is.null(p)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Plot could not be generated.")
        return()
      }

      print(p)
    }
    
  )
)
