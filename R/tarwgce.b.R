TARWGCEClass <- R6::R6Class(
  "TARWGCEClass",
  inherit = TARWGCEBase,
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
      
      summaryRows <- c(
        formula = "Formula",
        n = "Complete cases",
        cvNfolds = "Cross-validation folds",
        errorMeasure = "Prediction-error measure",
        errorMeasureWhich = "Selection rule",
        seed = "Seed",
        method = "Method",
        supportMethod = "Support method",
        signalSupportSpecification = "Signal support specification",
        ridgeLambdaRange = "Ridge lambda range",
        ridgeLambdaN = "Number of ridge lambda values",
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
      
      summaryTable <- self$results$modelSummary
      for (key in names(summaryRows)) {
        summaryTable$addRow(
          rowKey = key,
          values = list(measure = summaryRows[[key]])
        )
      }
      
      weightValues <- suppressWarnings(as.numeric(
        trimws(strsplit(self$options$weight, ",", fixed = TRUE)[[1L]])
      ))
      weightValues <- unique(weightValues[is.finite(weightValues)])
      
      self$results$plotCV$setSize(
        400,
        160 + 140 * max(1L, length(weightValues))
      )
      
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

      self$results$plotRidge$setVisible(
        isTRUE(self$options$plotRidge))
      
      self$results$plot1$setVisible(
        isTRUE(self$options$plot1))
      
      self$results$plot2$setVisible(
        isTRUE(self$options$plot2) && supportGridAvailable)
      
      self$results$plot3$setVisible(
        isTRUE(self$options$plot3) && supportGridAvailable)
      
      self$results$plot4$setVisible(
        isTRUE(self$options$plot4) && supportGridAvailable)
      
      self$results$plot5$setVisible(
        isTRUE(self$options$plot5) && supportGridAvailable)
      
      self$results$plot6$setVisible(
        isTRUE(self$options$plot6) && reestimationPlotsAvailable)
      
      self$results$plot7$setVisible(
        isTRUE(self$options$plot7) && reestimationPlotsAvailable)
      
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
      
      formStr <- jmvcore::constructFormula(
        dep = dep,
        terms = c(covs, factors)
      )
      
      if (!isTRUE(self$options$intercept))
        formStr <- paste(formStr, "- 1")
      
      form <- stats::as.formula(formStr)
      
      if (needsTrueCoef) {
        
        trueCoefText <- trimws(self$options$trueCoef)
        
        if (trueCoefText == "") {
          jmvcore::reject(
            "The precision-error plots require a vector of true coefficients.",
            code = "true_coef_required"
          )
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
      
      if (
        self$options$ridgeLambdaMin >= self$options$ridgeLambdaMax
      ) {
        jmvcore::reject("The minimum ridge lambda must be smaller than the maximum ridge lambda.", code = "invalid_lambda_range")
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
        isTRUE(self$options$plotCV) && cvPlotAvailable)
      
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
      useCachedFit <- is.list(fitCache) &&
        !is.null(fitCache$fit) && !is.null(fitCache$ci)

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
          jmvcore::reject("The signal support value must be a positive finite number.",
            code = "invalid_signal_support")
        }
      }
      
      if (useCachedFit) {
        fit <- fitCache$fit
        warningMsg <- fitCache$warningMsg
        fitElapsed <- fitCache$fitElapsed
      } else {
        fit <- tryCatch(
          withCallingHandlers(
          GCEstim::cv.lmgce(
            formula = form,
            data = data,
            support.method = "ridge",
            support.method.ridge.lambda.min = self$options$ridgeLambdaMin,
            support.method.ridge.lambda.max = self$options$ridgeLambdaMax,
            support.method.ridge.lambda.n = self$options$ridgeLambdaN,
            support.method.ridge.symm = 
              (self$options$ridgeSignalSupport == "symm"),
            support.method.ridge.maxresid = 
              (self$options$noiseSupportMethod == "maxres"),
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
            OLS = self$options$OLS,
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
      
      elapsed <- if (useCachedFit) fitElapsed else as.numeric(
        difftime(Sys.time(), start_time, units = "secs")
      )
      fitElapsed <- elapsed
      
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
        ci <- fitCache$ci
      } else {
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
        

        self$results$fitCache$setState(list(
          fit = fit,
          ci = ci,
          warningMsg = fitWarningMsg,
          fitElapsed = fitElapsed
        ))
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
        measure = "Method",
        value = methodLabel
      ))
      
      summaryTable$setRow(rowKey = "supportMethod", values = list(
        measure = "Support method",
        value = "Ridge"
      ))
      
      summaryTable$setRow(rowKey = "signalSupportSpecification", values = list(
        measure = "Signal support specification",
        value = if (self$options$ridgeSignalSupport == "symm") {
          "Symmetric"
        } else {
          "Asymmetric"
        }
      ))
      
      summaryTable$setRow(rowKey = "ridgeLambdaRange", values = list(
        measure = "Ridge lambda range",
        value = paste0(
          "[",
          self$options$ridgeLambdaMin,
          ", ",
          self$options$ridgeLambdaMax,
          "]"
        )
      ))
      
      summaryTable$setRow(rowKey = "ridgeLambdaN", values = list(
        measure = "Number of ridge lambda values",
        value = as.character(self$options$ridgeLambdaN)
      ))
      
      noiseResponse <- as.numeric(data[[dep]])
      
      ridgeFit <- if (useCachedFit) fitCache$ridgeFit else NULL
      if (self$options$noiseSupportMethod == "sigma3") {
        noiseLimit <- 3 * stats::sd(noiseResponse)
      } else if (useCachedFit && !is.null(fitCache$noiseLimit)) {
        noiseLimit <- fitCache$noiseLimit
      } else {
        ridgeFit <- tryCatch(
          GCEstim::ridgetrace(
            formula = form,
            data = data,
            lambda.min = self$options$ridgeLambdaMin,
            lambda.max = self$options$ridgeLambdaMax,
            lambda.n = self$options$ridgeLambdaN,
            errormeasure = self$options$errorMeasure,
            cv = FALSE,
            seed = self$options$seed
          ),
          error = function(e) e
        )
        if (inherits(ridgeFit, "error")) {
          jmvcore::reject("Error computing noise support: {message}",
            code = "noise_support_failed", message = conditionMessage(ridgeFit))
        }
        noiseLimit <- ridgeFit$max.abs.residual
      }

      if (!useCachedFit || is.null(fitCache$noiseLimit)) {
        self$results$fitCache$setState(list(
          fit = fit,
          ci = ci,
          warningMsg = fitWarningMsg,
          fitElapsed = fitElapsed,
          noiseLimit = noiseLimit,
          ridgeFit = ridgeFit
        ))
      }
      
      noiseInterval <- paste0(
        "[",
        format(signif(-noiseLimit, 6), trim = TRUE), ", ",
        format(signif(noiseLimit, 6), trim = TRUE), "]"
      )
      
      summaryTable$setRow(rowKey = "noiseSupportMethod", values = list(
        measure = "Noise support specification",
        value = paste0(
          if (self$options$noiseSupportMethod == "sigma3")
            "3 sigma: "
          else
            "Maximum absolute ridge residuals: ",
          noiseInterval
        )
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
          "Asymptotic (z)"
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
      
      cvTable <- self$results$cvResults
      
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
      
      plotWarnings <- character()
      
      tryPlot <- function(label, expr) {
        tryCatch(
          expr,
          error = function(e) {
            plotWarnings <<- c(
              plotWarnings,
              paste0(label, ": ", conditionMessage(e))
            )
            NULL
          }
        )
      }
      
      plotToPNG <- function(p, height = 900) {
        path <- tempfile(fileext = ".png")
        on.exit(unlink(path), add = TRUE)
        
        grDevices::png(path, width = 1200, height = height, res = 150)
        tryCatch(print(p), finally = grDevices::dev.off())
        
        readBin(path, what = "raw", n = file.info(path)$size)
      }
      
      # Retained PNG states avoid regenerating unchanged plots.
      needsPlot <- function(result) {
        state <- result$state
        is.null(state) || is.null(state$png)
      }

      plots <- list()

      if (self$options$plot1 && needsPlot(self$results$plot1))
        plots$plot1 <- tryPlot(
          "Estimates and confidence intervals",
          plot(
            bestFit,
            which = 1,
            ci.level = self$options$bootConfLevel,
            ci.method = ifelse(
              self$options$bootstrap,
              self$options$bootCIMethod,
              "z"
            ),
            OLS = self$options$OLS
          )
        )

      if (self$options$plot2 && supportGridAvailable && needsPlot(self$results$plot2))
        plots$plot2 <- tryPlot(
          "Supports vs prediction error",
          plot(
            bestFit,
            which = 2,
            OLS = self$options$OLS,
            NormEnt = self$options$NormEnt
          )
        )
      
      if (self$options$plot3 && supportGridAvailable && needsPlot(self$results$plot3))
        plots$plot3 <- tryPlot(
          "Supports vs estimates",
          plot(
            bestFit,
            which = 3,
            OLS = self$options$OLS
          )
        )
      
      if (self$options$plot4 && supportGridAvailable && needsPlot(self$results$plot4))
        plots$plot4 <- tryPlot(
          "Supports vs normalized entropy",
          plot(
            bestFit,
            which = 4,
            OLS = self$options$OLS
          )
        )
      
      if (self$options$plot5 && supportGridAvailable && needsPlot(self$results$plot5))
        plots$plot5 <- tryPlot(
          "Supports vs precision error",
          plot(
            bestFit,
            which = 5,
            coef = trueCoef,
            OLS = self$options$OLS,
            NormEnt = self$options$NormEnt
          )
        )
      
      if (self$options$plot6 && reestimationPlotsAvailable && needsPlot(self$results$plot6))
        plots$plot6 <- tryPlot(
          "GCE re-estimation vs prediction error",
          plot(
            bestFit,
            which = 6,
            OLS = self$options$OLS
          )
        )
      
      if (self$options$plot7 && reestimationPlotsAvailable && needsPlot(self$results$plot7))
        plots$plot7 <- tryPlot(
          "GCE re-estimation vs precision error",
          plot(
            bestFit,
            which = 7,
            coef = trueCoef,
            OLS = self$options$OLS
          )
        )
      
      if (isTRUE(self$options$plotCV) && cvPlotAvailable && needsPlot(self$results$plotCV))
        plots$plotCV <- tryPlot(
          "Cross-validation error by support combination",
          plot(
            fit,
            which = 1,
            ncol = 1,
            scales = "free"
          )
        )
      
      if (self$options$plotRidge && needsPlot(self$results$plotRidge)) {
        plots$plotRidge <- tryPlot(
          "Ridge trace",
          {
            if (is.null(ridgeFit))
              ridgeFit <- GCEstim::ridgetrace(
                formula = form,
                data = data,
                lambda.min = self$options$ridgeLambdaMin,
                lambda.max = self$options$ridgeLambdaMax,
                lambda.n = self$options$ridgeLambdaN,
                errormeasure = self$options$errorMeasure,
                cv = FALSE,
                seed = self$options$seed
              )
            
            plot(ridgeFit, log = TRUE, range = TRUE)
          }
        )
      }

      firstPlot <- function(x) {
        if (inherits(x, "ggplot"))
          return(x)
        
        if (is.null(x) || length(x) == 0L)
          return(NULL)
        
        if (inherits(x[[1L]], "ggplot"))
          return(x[[1L]])
        
        NULL
      }
      
      setPlotPNG <- function(result, x, label, height = 900) {
        p <- firstPlot(x)
        
        bytes <- if (is.null(p)) {
          NULL
        } else {
          tryPlot(label, plotToPNG(p, height = height))
        }
        
        result$setState(list(png = bytes))
      }
      
      if (self$options$plot1 && needsPlot(self$results$plot1)) {
        p <- plots$plot1[["p1"]]
        
        pngBytes <- if (is.null(p)) {
          NULL
        } else {
          tryPlot("Estimates and confidence intervals", plotToPNG(p))
        }
        
        self$results$plot1$setState(list(png = pngBytes))
      }

      if (self$options$plot2 && supportGridAvailable && needsPlot(self$results$plot2))
        setPlotPNG(
          self$results$plot2, plots$plot2,
          "Supports vs prediction error"
        )
      
      if (self$options$plot3 && supportGridAvailable && needsPlot(self$results$plot3))
        setPlotPNG(
          self$results$plot3, plots$plot3,
          "Supports vs estimates"
        )
      
      if (self$options$plot4 && supportGridAvailable && needsPlot(self$results$plot4))
        setPlotPNG(
          self$results$plot4, plots$plot4,
          "Supports vs normalized entropy"
        )
      
      if (self$options$plot5 && supportGridAvailable && needsPlot(self$results$plot5))
        setPlotPNG(
          self$results$plot5, plots$plot5,
          "Supports vs precision error"
        )
      
      if (self$options$plot6 && reestimationPlotsAvailable && needsPlot(self$results$plot6))
        setPlotPNG(
          self$results$plot6, plots$plot6,
          "GCE re-estimation vs prediction error"
        )
      
      if (self$options$plot7 && reestimationPlotsAvailable && needsPlot(self$results$plot7))
        setPlotPNG(
          self$results$plot7, plots$plot7,
          "GCE re-estimation vs precision error"
        )
      
      if (isTRUE(self$options$plotCV) && cvPlotAvailable && needsPlot(self$results$plotCV))
        setPlotPNG(
          self$results$plotCV, plots$plotCV,
          "Cross-validation error by support combination",
          height = round(3 * (160 + 140 * length(weightValues)))
        )
      
      if (self$options$plotRidge && needsPlot(self$results$plotRidge))
        setPlotPNG(
          self$results$plotRidge, plots$plotRidge,
          "Ridge trace"
        )
      
      ## Complete
      
      warningMessages <- plotWarnings
      
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
    
    .plot1 = function(image, ggtheme, theme) {
      if (is.null(image$state) || is.null(image$state$png))
        return(FALSE)
      
      grid::grid.newpage()
      grid::grid.raster(
        png::readPNG(image$state$png, native = TRUE),
        width = grid::unit(1, "npc"),
        height = grid::unit(1, "npc"),
        interpolate = TRUE
      )
      
      TRUE
    },
    .plot2 = function(image, ggtheme, theme) private$.printPNG(image),
    .plot3 = function(image, ggtheme, theme) private$.printPNG(image),
    .plot4 = function(image, ggtheme, theme) private$.printPNG(image),
    .plot5 = function(image, ggtheme, theme) private$.printPNG(image),
    .plot6 = function(image, ggtheme, theme) private$.printPNG(image),
    .plot7 = function(image, ggtheme, theme) private$.printPNG(image),
    .plotCV = function(image, ggtheme, theme) private$.printPNG(image),
    .plotRidge = function(image, ggtheme, theme) private$.printPNG(image),

    
    .printPNG = function(image) {
      if (is.null(image$state) || is.null(image$state$png))
        return(FALSE)
      
      grid::grid.newpage()
      grid::grid.raster(
        png::readPNG(image$state$png, native = TRUE),
        width = grid::unit(1, "npc"),
        height = grid::unit(1, "npc"),
        interpolate = TRUE
      )
      
      TRUE
    }
  )
)
