LMGCEClass <- R6::R6Class(
  "LMGCEClass",
  inherit = LMGCEBase,
  private = list(
    
    .run = function() {
      
      # .libPaths(c(
      #   "C:/Users/Utilizador/AppData/Local/R/win-library/4.5",
      #   .libPaths()
      # ))
      
      .libPaths(c(
            "C:/Users/jorge/AppData/Local/R/win-library/4.5",
            .libPaths()
       ))
        
      dep  <- unlist(self$options$dep)
      covs    <- unlist(self$options$covs)
      factors <- unlist(self$options$factors)
      
      if (self$options$run == 0) {
        self$results$text$setContent(
          "<p>Select variables and options, then press <b>Run model</b>.</p>"
        )
        return()
      }
      
      if (length(dep) == 0 || (length(covs) == 0 && length(factors) == 0)) {
        self$results$text$setContent(
          "<p>Select one dependent variable and at least one covariate or factor.</p>"
        )
        return()
      }
      
      selectedPlots <- c(
        self$options$plot1,
        self$options$plot2,
        self$options$plot3,
        self$options$plot4,
        self$options$plot5,
        self$options$plot6,
        self$options$plot7,
        self$options$plotRidge &&
          self$options$supportMethod == "ridge"
      )
      
      if (sum(selectedPlots) > 3) {
        self$results$text$setContent(
          "<p style='color:red;'><b>Error:</b> Please select a maximum of 3 plots.</p>"
        )
        return()
      }
      
      needsTrueCoef <- self$options$plot5 || self$options$plot7
      
      trueCoef <- NULL
      
      vars <- c(dep, covs, factors)
      
      data <- self$data[, vars, drop = FALSE]
      data <- data[stats::complete.cases(data), , drop = FALSE]
      
      if (nrow(data) < 3) {
        self$results$text$setContent(
          "<p>Not enough complete cases to estimate the model.</p>"
        )
        return()
      }
      
      for (f in factors) {
        data[[f]] <- as.factor(data[[f]])
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
          self$results$text$setContent(
            "<p style='color:red;'><b>Error:</b> Plots 5 and 7 require a vector of true coefficients.</p>"
          )
          return()
        }
        
        trueCoef <- suppressWarnings(
          as.numeric(strsplit(trueCoefText, ",")[[1]])
        )
        
        if (any(is.na(trueCoef))) {
          self$results$text$setContent(
            "<p style='color:red;'><b>Error:</b> True coefficients must be numeric and separated by commas.</p>"
          )
          return()
        }
        
        X <- stats::model.matrix(form, data)
        
        expectedCoef <- ncol(X)
        
        if (length(trueCoef) != expectedCoef) {
          self$results$text$setContent(
            paste0(
              "<p style='color:red;'><b>Error:</b> The true coefficient vector must contain ",
              expectedCoef,
              " values (including the intercept).</p>"
            )
          )
          return()
        }
      }
      
      if (
        self$options$supportMethod == "ridge" &&
        self$options$ridgeLambdaMin >= self$options$ridgeLambdaMax
      ) {
        self$results$text$setContent(
          "<p style='color:red;'><b>Error:</b> The minimum ridge lambda must be smaller than the maximum ridge lambda.</p>"
        )
        return()
      }
      
      if (self$options$supportSignalVectorMin >= self$options$supportSignalVectorMax) {
        self$results$text$setContent(
          "<p style='color:red;'><b>Error:</b> The minimum signal support range must be smaller than the maximum.</p>"
        )
        return()
      }
      
      if (self$options$cvNfolds > nrow(data)) {
        self$results$text$setContent(
          paste0(
            "<p style='color:red;'><b>Error:</b> The number of CV folds cannot exceed the number of complete cases (",
            nrow(data),
            ").</p>"
          )
        )
        return()
      }
      
      if (self$options$M < 3 || self$options$M > 99 || self$options$M %% 2 == 0) {
        self$results$text$setContent(
          "<p style='color:red;'><b>Error:</b> Signal support points must be an odd integer between 3 and 99.</p>"
        )
        return()
      }
      
      if (self$options$J < 3 || self$options$J > 99 || self$options$J %% 2 == 0) {
        self$results$text$setContent(
          "<p style='color:red;'><b>Error:</b> Noise support points must be an odd integer between 3 and 99.</p>"
        )
        return()
      }
      
      bootB <- 0
      bootMethod <- "residuals"
      
      if (self$options$bootstrap) {
        
        bootB <- self$options$bootB
        bootMethod <- self$options$bootMethod
        
        if (bootB < 10) {
          self$results$text$setContent(
            "<p style='color:red;'><b>Error:</b> Number of bootstrap samples must be at least 10.</p>"
          )
          return()
        }
      }
      
      self$results$text$setContent("<p>Starting fit...</p>")
      
      start_time <- Sys.time()
      
      warningMsg <- NULL
      
      fit <- tryCatch(
        withCallingHandlers(
        GCEstim::lmgce(
          formula = form,
          data = data,
          support.method = self$options$supportMethod,
          support.method.ridge.lambda.min = self$options$ridgeLambdaMin,
          support.method.ridge.lambda.max = self$options$ridgeLambdaMax,
          support.method.ridge.lambda.n = self$options$ridgeLambdaN,
          support.method.ridge.symm =
            (self$options$supportMethod != "ridge" ||
            self$options$ridgeSignalSupport == "symm"),
          support.method.ridge.maxresid =
            (self$options$supportMethod == "ridge" &&
            self$options$noiseSupportMethod == "maxres"),
          support.signal.vector.n = self$options$supportSignalVectorN,
          support.signal.vector.min = self$options$supportSignalVectorMin,
          support.signal.vector.max = self$options$supportSignalVectorMax,
          support.signal.points = self$options$M,
          support.noise.points = self$options$J,
          cv.nfolds = self$options$cvNfolds,
          errormeasure = self$options$errorMeasure,
          errormeasure.which = self$options$errorMeasureWhich,
          weight = self$options$weight,
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
        self$results$text$setContent(
          paste0("<p><b>Error:</b> ", fit$message, "</p>")
        )
        return()
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
      
      ## Confidence Interval
      
      ci <- NULL
      
      if (self$options$bootstrap) {
        
        ci <- tryCatch(
          stats::confint(
            fit,
            level = self$options$bootConfLevel,
            method = self$options$bootCIMethod
          ),
          error = function(e) e
        )
        
      } else {
        
        ci <- tryCatch(
          stats::confint(
            fit,
            level = self$options$bootConfLevel
          ),
          error = function(e) e
        )
      }
      
      if (inherits(ci, "error")) {
        self$results$text$setContent(
          paste0(
            "<p style='color:red;'><b>Error in confidence intervals:</b> ",
            ci$message,
            "</p>"
          )
        )
        return()
      }
      
      self$results$text$setContent("<p>CI completed. Filling tables...</p>")
      
      ## Summary table
      
      summ <- tryCatch(
        summary(fit),
        error = function(e) e
      )
      
      if (inherits(summ, "error")) {
        self$results$text$setContent(
          paste0("<p><b>Error in summary:</b> ", summ$message, "</p>")
        )
        return()
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
        value = self$options$supportMethod
      ))
      
      summaryTable$addRow(rowKey = "signalSupportSpecification", values = list(
        measure = "Signal support specification",
        value = if (
          self$options$supportMethod == "standardized" ||
          self$options$ridgeSignalSupport == "symm"
        ) {
          "Symmetric"
        } else {
          "Asymmetric"
        }
      ))
      
      if (self$options$supportMethod == "ridge") {
        
        summaryTable$addRow(rowKey = "ridgeLambdaRange", values = list(
          measure = "Ridge lambda range",
          value = paste0(
            "[",
            self$options$ridgeLambdaMin,
            ", ",
            self$options$ridgeLambdaMax,
            "]"
          )
        ))
        
        summaryTable$addRow(rowKey = "ridgeLambdaN", values = list(
          measure = "Number of ridge lambda values",
          value = as.character(self$options$ridgeLambdaN)
        ))
      }
      
      summaryTable$addRow(rowKey = "noiseSupportMethod", values = list(
        measure = "Noise support specification",
        value = if (
          self$options$supportMethod == "standardized" ||
          self$options$noiseSupportMethod == "sigma3"
        ) {
          "3 sigma"
        } else {
          "Maximum absolute ridge residuals"
        }
      ))
      
      summaryTable$addRow(rowKey = "supportSignalVectorN", values = list(
        measure = "Number of support spaces",
        value = as.character(self$options$supportSignalVectorN)
      ))
      
      summaryTable$addRow(rowKey = "M", values = list(
        measure = "Signal support points",
        value = as.character(self$options$M)
      ))
      
      summaryTable$addRow(rowKey = "J", values = list(
        measure = "Noise support points",
        value = as.character(self$options$J)
      ))
      
      summaryTable$addRow(rowKey = "weight", values = list(
        measure = "Noise weight",
        value = as.character(self$options$weight)
      ))
      
      summaryTable$addRow(rowKey = "caseGLM", values = list(
        measure = "GLM case",
        value = self$options$caseGLM
      ))
      
      summaryTable$addRow(rowKey = "twostepsN", values = list(
        measure = "Post-GCE reestimations",
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
            ciLowerZ = if (!hasCI) ci[term, 1] else NA,
            ciUpperZ = if (!hasCI) ci[term, 2] else NA,
            ciLower = if (hasCI) ci[term, 1] else NA,
            ciUpper = if (hasCI) ci[term, 2] else NA,
            se = coefs[i, "Std. Deviation"],
            z = coefs[i, "z value"],
            p = coefs[i, "Pr(>|t|)"]
          )
        )
      }
      
      self$results$text$setContent("<p>Tables completed. Plotting...</p>")
      
      plots <- list()

      if (self$options$plot1)
        plots$plot1 <- tryCatch(plot(fit,
                                     which = 1,
                                     ci.level = self$options$bootConfLevel,
                                     ci.method = ifelse(self$options$bootstrap,
                                                        self$options$bootCIMethod,
                                                        "z"),
                                     OLS = self$options$OLS),
                                error = function(e) NULL)

      if (self$options$plot2)
        plots$plot2 <- tryCatch(plot(fit,
                                     which = 2,
                                     OLS = self$options$OLS,
                                     NormEnt = self$options$NormEnt),
                                error = function(e) NULL)

      if (self$options$plot3)
        plots$plot3 <- tryCatch(plot(fit,
                                     which = 3,
                                     OLS = self$options$OLS),
                                error = function(e) NULL)

      if (self$options$plot4)
        plots$plot4 <- tryCatch(plot(fit,
                                     which = 4,
                                     OLS = self$options$OLS),
                                error = function(e) NULL)

      if (self$options$plot5)
        plots$plot5 <- tryCatch(plot(fit,
                                     which = 5,
                                     coef = trueCoef,
                                     OLS = self$options$OLS,
                                     NormEnt = self$options$NormEnt),
                                error = function(e) NULL)

      if (self$options$plot6)
        plots$plot6 <- tryCatch(plot(fit,
                                     which = 6,
                                     OLS = self$options$OLS),
                                error = function(e) NULL)

      if (self$options$plot7)
        plots$plot7 <- tryCatch(plot(fit,
                                     which = 7,
                                     coef = trueCoef,
                                     OLS = self$options$OLS),
                                error = function(e) NULL)
      
      if (
        self$options$plotRidge &&
        self$options$supportMethod == "ridge"
      ) {
        plots$plotRidge <- tryCatch(
          {
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
            
            plot(
              ridgeFit,
              log = TRUE,
              range = TRUE
            )
          },
          error = function(e) NULL
        )
      }

      if (self$options$plot1)
        self$results$plot1$setState(list(p = plots$plot1))

      if (self$options$plot2)
        self$results$plot2$setState(list(p = plots$plot2))

      if (self$options$plot3)
        self$results$plot3$setState(list(p = plots$plot3))

      if (self$options$plot4)
        self$results$plot4$setState(list(p = plots$plot4))

      if (self$options$plot5)
        self$results$plot5$setState(list(p = plots$plot5))

      if (self$options$plot6)
        self$results$plot6$setState(list(p = plots$plot6))

      if (self$options$plot7)
        self$results$plot7$setState(list(p = plots$plot7))
      
      if (
        self$options$plotRidge &&
        self$options$supportMethod == "ridge"
      )
        self$results$plotRidge$setState(
          list(p = plots$plotRidge)
        )
      
      ## Complete
      
      if (!is.null(warningMsg)) {
        
        self$results$text$setContent(
          paste0(
            "<p><b>GCE linear model estimated successfully.</b></p>",
            "<p><b>Computation time:</b> ", time_txt, "</p>",
            "<p style='color:orange;'><b>Warning:</b> ",
            warningMsg,
            "</p>"
          )
        )
        
      } else {
        
        self$results$text$setContent(
          paste0(
            "<p>GCE linear model estimated successfully.</p>",
            "<p><b>Computation time:</b> ", time_txt, "</p>"
          )
        )
        
      }
      
      ## References
      
      self$results$reference$setContent(
        paste0(
          "<p>Cabral, J. (2026). <i>GCEstim: Regression Coefficients Estimation ",
          "Using the Generalized Cross Entropy</i>. R package version ",
          as.character(utils::packageVersion("GCEstim")),
          ". https://CRAN.R-project.org/package=GCEstim</p>"
        )
      )
    },
    
    .plot1 = function(image, ggtheme, theme) private$.printPlot(image),
    .plot2 = function(image, ggtheme, theme) private$.printPlot(image),
    .plot3 = function(image, ggtheme, theme) private$.printPlot(image),
    .plot4 = function(image, ggtheme, theme) private$.printPlot(image),
    .plot5 = function(image, ggtheme, theme) private$.printPlot(image),
    .plot6 = function(image, ggtheme, theme) private$.printPlot(image),
    .plot7 = function(image, ggtheme, theme) private$.printPlot(image),
    .plotRidge = function(image, ggtheme, theme) private$.printPlot(image),

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