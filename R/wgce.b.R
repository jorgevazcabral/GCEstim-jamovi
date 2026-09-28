WGCEClass <- R6::R6Class(
  "WGCEClass",
  inherit = WGCEBase,
  private = list(

    .run = function() {
      
      # .libPaths(c(
      #   "C:/Users/jorge/AppData/Local/R/win-library/4.5",
      #   .libPaths()
      # ))

      escapeHTML <- function(x) {
        x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
        x <- gsub("<", "&lt;", x, fixed = TRUE)
        x <- gsub(">", "&gt;", x, fixed = TRUE)
        x
      }

      showError <- function(message) {
        self$results$text$setContent(paste0(
          "<p style='color:red;'><b>Error:</b> ",
          escapeHTML(message), "</p>"
        ))
      }

      parseLimits <- function(input, label) {
        input <- trimws(input)
        if (!nzchar(input))
          stop(paste0(label, " cannot be empty."), call. = FALSE)
        parts <- trimws(strsplit(input, ",", fixed = TRUE)[[1]])
        if (any(!nzchar(parts)))
          stop(paste0(label, " must contain numbers separated by commas."),
               call. = FALSE)
        values <- suppressWarnings(as.numeric(parts))
        if (any(!is.finite(values)))
          stop(paste0(label, " must contain only finite numbers."),
               call. = FALSE)
        values
      }

      if (self$options$run == 0) {
        self$results$text$setContent(
          "<p>Select variables and support limits, then press <b>Run model</b>.</p>"
        )
        return()
      }

      dep <- unlist(self$options$dep)
      covs <- unlist(self$options$covs)
      factors <- unlist(self$options$factors)

      if (length(dep) != 1L ||
          (length(covs) == 0L && length(factors) == 0L)) {
        showError("Select one dependent variable and at least one covariate or factor.")
        return()
      }

      vars <- unique(c(dep, covs, factors))
      data <- self$data[, vars, drop = FALSE]
      data <- data[stats::complete.cases(data), , drop = FALSE]

      if (nrow(data) < 3L) {
        showError("Not enough complete cases to estimate the model.")
        return()
      }

      for (f in factors)
        data[[f]] <- droplevels(as.factor(data[[f]]))

      form <- stats::reformulate(
        termlabels = c(covs, factors),
        response = dep,
        intercept = self$options$intercept
      )

      X <- tryCatch(stats::model.matrix(form, data), error = function(e) e)
      if (inherits(X, "error")) {
        showError(conditionMessage(X))
        return()
      }
      coefNames <- colnames(X)
      k <- length(coefNames)

      if (k == 0L) {
        showError("The model has no coefficients.")
        return()
      }

      if (self$options$M %% 2L != 1L ||
          self$options$J %% 2L != 1L) {
        showError("Signal and noise support points must each be odd integers.")
        return()
      }

      limits <- tryCatch(list(
        lower = parseLimits(self$options$supportSignalLowers,
                            "Lower limits"),
        upper = parseLimits(self$options$supportSignalUppers,
                            "Upper limits")
      ), error = function(e) e)

      if (inherits(limits, "error")) {
        showError(conditionMessage(limits))
        return()
      }

      if (!(length(limits$lower) %in% c(1L, k)) ||
          !(length(limits$upper) %in% c(1L, k))) {
        showError(paste0(
          "Enter either one value or ", k, " values in each field. ",
          "Coefficient order: ", paste(coefNames, collapse = ", "), "."
        ))
        return()
      }

      lower <- rep(limits$lower, length.out = k)
      upper <- rep(limits$upper, length.out = k)

      if (any(!is.finite(lower)) || any(!is.finite(upper)) ||
          any(lower >= upper)) {
        showError("Each signal support must have a finite lower limit smaller than its upper limit.")
        return()
      }

      # A length-two vector shares its interval; a k-by-two matrix assigns
      # one interval to each column of the model matrix, including intercept.
      supportSignal <- if (length(limits$lower) == 1L &&
                           length(limits$upper) == 1L) {
        c(lower[1], upper[1])
      } else {
        matrix(c(lower, upper), nrow = k, ncol = 2L,
               dimnames = list(coefNames, c("LL", "UL")))
      }

      supportTable <- self$results$signalSupports
      for (i in seq_len(k)) {
        supportTable$addRow(rowKey = as.character(i), values = list(
          term = coefNames[i], lower = lower[i], upper = upper[i]
        ))
      }

      supportNoise <- NULL
      if (!self$options$supportNoise3sig) {
        if (self$options$supportNoiseMin >= self$options$supportNoiseMax) {
          showError("Minimum noise support must be smaller than maximum noise support.")
          return()
        }
        supportNoise <- c(self$options$supportNoiseMin,
                          self$options$supportNoiseMax)
      }

      bootB <- if (self$options$bootstrap) self$options$bootB else 0L

      self$results$text$setContent("<p>Estimating W-GCE model...</p>")
      startTime <- Sys.time()
      warnings <- character()

      fit <- tryCatch(
        withCallingHandlers(
          GCEstim::lmgce(
            formula = form,
            data = data,
            cv = FALSE,
            errormeasure.which = "min",
            support.method = "standardized",
            support.signal = supportSignal,
            support.signal.points = self$options$M,
            support.noise = supportNoise,
            support.noise.points = self$options$J,
            weight = self$options$weight,
            twosteps.n = 0L,
            method = self$options$method,
            caseGLM = self$options$caseGLM,
            boot.B = bootB,
            boot.method = self$options$bootMethod,
            seed = self$options$seed,
            OLS = self$options$OLS
          ),
          warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) e
      )

      if (inherits(fit, "error")) {
        showError(conditionMessage(fit))
        return()
      }

      if (self$options$saveFitted &&
          self$results$saveFitted$isNotFilled()) {
        self$results$saveFitted$setRowNums(rownames(data))
        self$results$saveFitted$setValues(as.numeric(fit$fitted.values))
      }

      if (self$options$saveResiduals &&
          self$results$saveResiduals$isNotFilled()) {
        self$results$saveResiduals$setRowNums(rownames(data))
        self$results$saveResiduals$setValues(as.numeric(fit$residuals))
      }

      if (self$options$saveW) {
        w <- as.matrix(fit$w)
        if (self$options$caseGLM != "D") {
          warnings <- c(warnings,
                        "Noise probabilities cannot be saved as observation columns for the moment-based GLM cases.")
        } else if (nrow(w) != nrow(data) || ncol(w) != self$options$J) {
          warnings <- c(warnings,
                        "Noise probabilities do not match the expected observation and support-point dimensions.")
        } else if (self$results$saveW$isNotFilled()) {
          self$results$saveW$setRowNums(rownames(data))
          for (j in seq_len(ncol(w))) {
            self$results$saveW$setValues(as.numeric(w[, j]), index = j)
          }
        }
      }

      elapsed <- as.numeric(difftime(Sys.time(), startTime,
                                     units = "secs"))
      timeText <- if (elapsed < 60) {
        paste0(round(elapsed, 3), " s")
      } else {
        paste0(floor(elapsed / 60), " min ",
               round(elapsed %% 60, 1), " s")
      }

      ci <- tryCatch(
        if (self$options$bootstrap) {
          stats::confint(fit, level = self$options$bootConfLevel,
                         method = self$options$bootCIMethod)
        } else {
          stats::confint(fit, level = self$options$bootConfLevel)
        },
        error = function(e) e
      )
      if (inherits(ci, "error")) {
        showError(paste0("Confidence intervals: ", conditionMessage(ci)))
        return()
      }

      summ <- tryCatch(summary(fit), error = function(e) e)
      if (inherits(summ, "error")) {
        showError(paste0("Model summary: ", conditionMessage(summ)))
        return()
      }

      summaryTable <- self$results$modelSummary
      addMeasure <- function(key, label, value) {
        summaryTable$addRow(rowKey = key,
                            values = list(measure = label,
                                          value = as.character(value)))
      }

      addMeasure("formula", "Formula", paste(deparse(form), collapse = " "))
      addMeasure("n", "Complete cases", nrow(data))
      addMeasure("supportMode", "Signal support", if (
        length(unique(lower)) == 1L && length(unique(upper)) == 1L
      ) paste0("[", lower[1], ", ", upper[1], "]") else
        "Different for each coefficient (see Signal Supports)")
      addMeasure("M", "Signal support points", self$options$M)
      addMeasure("J", "Noise support points", self$options$J)
      addMeasure("noise", "Noise support", if (
        self$options$supportNoise3sig
      ) {
        noiseResponse <- as.numeric(data[[dep]])
        if (self$options$caseGLM %in% c("M", "NM")) {
          noiseResponse <- as.numeric(crossprod(X, noiseResponse))
          if (self$options$caseGLM == "NM")
            noiseResponse <- noiseResponse / nrow(X)
        }
        noiseLimit <- 3 * stats::sd(noiseResponse)
        paste0(
          "Automatic (3 sigma): [",
          format(signif(-noiseLimit, 6), trim = TRUE), ", ",
          format(signif(noiseLimit, 6), trim = TRUE), "]"
        )
      } else paste0(
        "[", self$options$supportNoiseMin, ", ",
        self$options$supportNoiseMax, "]"
      ))
      addMeasure("weight", "Noise weight", self$options$weight)
      addMeasure("method", "Optimization method", self$options$method)
      addMeasure("case", "GLM case", self$options$caseGLM)
      addMeasure("boot", "Bootstrap", if (
        self$options$bootstrap
      ) "Yes" else "No")
      if (self$options$bootstrap) {
        addMeasure("bootB", "Bootstrap samples", self$options$bootB)
        addMeasure("bootMethod", "Bootstrap method", self$options$bootMethod)
      }
      addMeasure("ciLevel", "Confidence level",
                 paste0(100 * self$options$bootConfLevel, "%"))

      coefs <- summ$coefficients
      coefTable <- self$results$coefficients
      hasCI <- isTRUE(self$options$bootstrap)
      for (i in seq_len(nrow(coefs))) {
        term <- rownames(coefs)[i]
        coefTable$addRow(rowKey = term, values = list(
          term = term,
          estimate = coefs[i, "Estimate"],
          ciLowerZ = if (!hasCI) ci[term, 1] else NA_real_,
          ciUpperZ = if (!hasCI) ci[term, 2] else NA_real_,
          ciLower = if (hasCI) ci[term, 1] else NA_real_,
          ciUpper = if (hasCI) ci[term, 2] else NA_real_,
          se = coefs[i, "Std. Deviation"],
          z = coefs[i, "z value"],
          p = coefs[i, "Pr(>|t|)"]
        ))
      }

      
      pMatrix <- as.matrix(fit$p)
      if (nrow(pMatrix) != k || ncol(pMatrix) != self$options$M) {
        warnings <- c(warnings,
                      "Signal probabilities have unexpected dimensions and could not be displayed.")
      } else {
        pTable <- self$results$signalProbabilities
        pKeys <- paste0("p", seq_len(ncol(pMatrix)))
        for (j in seq_len(ncol(pMatrix))) {
          pTable$addColumn(
            name = pKeys[j], title = paste0("p_", j),
            type = "number", format = "zto"
          )
        }
        pTerms <- rownames(pMatrix)
        if (is.null(pTerms)) pTerms <- coefNames
        for (i in seq_len(nrow(pMatrix))) {
          values <- c(
            list(term = pTerms[i]),
            stats::setNames(as.list(as.numeric(pMatrix[i, ])), pKeys)
          )
          pTable$addRow(rowKey = as.character(i), values = values)
        }
      }

      if (self$options$plot1) {
        plot1 <- tryCatch(
          plot(fit, which = 1,
               ci.level = self$options$bootConfLevel,
               ci.method = if (hasCI) self$options$bootCIMethod else "z",
               OLS = self$options$OLS)[["p1"]],
          error = function(e) {
            warnings <<- c(warnings, paste0("Plot: ", conditionMessage(e)))
            NULL
          }
        )
        self$results$plot1$setState(list(p = plot1))
      }

      warningText <- if (length(warnings)) paste0(
        "<p style='color:orange;'><b>Warning:</b> ",
        escapeHTML(paste(unique(warnings), collapse = "; ")), "</p>"
      ) else ""
      self$results$text$setContent(paste0(
        "<p><b>W-GCE linear model estimated successfully.</b></p>",
        "<p><b>Computation time:</b> ", timeText, "</p>", warningText
      ))

      self$results$reference$setContent(paste0(
        "<p>Cabral, J. (2026). <i>GCEstim: Regression Coefficients ",
        "Estimation Using the Generalized Cross Entropy</i>. ",
        "R package version ",
        as.character(utils::packageVersion("GCEstim")),
        ". https://CRAN.R-project.org/package=GCEstim</p>"
      ))
    },

    .plot1 = function(image, ggtheme, theme) {
      p <- image$state$p
      if (is.null(p)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Plot could not be generated.")
      } else {
        print(p)
      }
    }
  )
)
