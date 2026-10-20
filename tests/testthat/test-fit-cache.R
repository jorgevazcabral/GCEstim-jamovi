cache_analysis <- function(kind, data, changes = list()) {
  ns <- asNamespace("GCEstimjmv")
  args <- modifyList(list(
    dep = "y", covs = c("X001", "X002"), run = TRUE,
    M = "3, 5", J = "3", weight = "0.5", cvNfolds = 2,
    supportSignalVectorN = 2, twostepsN = 1, bootstrap = FALSE,
    plot1 = FALSE, plot2 = FALSE, plotCV = FALSE, seed = 230676,
    trueCoef = "1, 1, 1"
  ), changes)
  if (kind == "TARWGCE") {
    args$ridgeLambdaN <- 100
    args$supportSignalVectorN <- 5
    args$cvNfolds <- 5
    args$J <- "3, 5"
    args$weight <- "0.3, 0.5"
  }
  options <- do.call(get(paste0(kind, "Options"), ns)$new, args)
  get(paste0(kind, "Class"), ns)$new(options = options, data = data)
}

cache_plain <- function(x) {
  if (inherits(x, c("lmgce", "ggplot")) || is.environment(x) ||
      is.function(x)) return(FALSE)
  if (is.list(x) && !all(vapply(x, cache_plain, logical(1)))) return(FALSE)
  attrs <- attributes(x)
  is.null(attrs) || all(vapply(attrs, cache_plain, logical(1)))
}

for (kind in c("TARWGCE", "TASWGCE")) local({
  kind <- kind
  testthat::test_that(paste(kind, "replays serialized cache and saved columns"), {
    data <- utils::read.csv(testthat::test_path("_fixtures", "dataThesis.csv"))
    data$X001[c(2, 9)] <- NA_real_
    original_cv <- GCEstim::cv.lmgce
    fullFit <- NULL
    testthat::local_mocked_bindings(cv.lmgce = function(...) {
      fullFit <<- original_cv(...)
      fullFit
    }, .package = "GCEstim")
    first <- cache_analysis(kind, data)
    first$run()
    state <- unserialize(serialize(first$results$fitCache$state, NULL))
    roundtrip <- cache_analysis(kind, data)
    roundtrip$results$fitCache$.__enclos_env__$super$fromProtoBuf(state_message(state),
                                          character(), character())
    testthat::expect_identical(roundtrip$results$fitCache$state, state)
    cached <- state
    testthat::expect_identical(cached$version, 4L)
    testthat::expect_true(cache_plain(cached))
    testthat::expect_true(all(c("plot1", "plot2", "plot3", "plot4", "plot6", "plotCV") %in%
                                names(cached$plots)))
    testthat::expect_true(all(vapply(cached$plots, function(x)
      is.data.frame(x$data) && nrow(x$data) > 0, logical(1))))
    testthat::expect_null(cached$payload)
    testthat::expect_true(cache_plain(cached$precisionBasis))
    testthat::expect_false(any(c("plot5", "plot7") %in% names(cached$plots)))
    clearWith <- first$results$fitCache$.__enclos_env__$private$.clearWith
    for (change in c("OLS", "NormEnt", "plot1", "plot5", "saveFitted",
                     "saveResiduals", "saveW", "cvTableRows", "trueCoef")) {
      testthat::expect_false(change %in% clearWith)
    }
    for (change in c("data", "M", "seed", "bootConfLevel")) {
      testthat::expect_true(change %in% clearWith)
    }

    # Old packed caches are ignored; the replacement is reused thereafter.
    calls <- 0L
    testthat::local_mocked_bindings(cv.lmgce = function(...) {
      calls <<- calls + 1L
      fullFit
    }, .package = "GCEstim")
    migrated <- cache_analysis(kind, data)
    legacy <- list(version = 3L, fit = state$fit, ci = state$ci,
                   coefficients = state$coefficients, plots = state$plots)
    migrated$results$fitCache$setState(list(payload =
      as.integer(memCompress(serialize(legacy, NULL), type = "gzip"))))
    migrated$run()
    migrated$run()
    testthat::expect_identical(calls, 1L)
    testthat::expect_identical(migrated$results$fitCache$state$version, 4L)

    # The reduced precision interface reproduces GCEstim's actual plots.
    for (which in c(5L, 7L)) {
      plotFit <- fullFit$best
      if (which == 5L) {
        sr <- plotFit$results$nocv$support.results
        plotFit$results$nocv$support.results <-
          sr[c(as.character(plotFit$support.ok), tail(names(sr), 3L))]
      }
      original <- plot(plotFit, which = which, coef = c(1, 1, 1),
                       OLS = TRUE, NormEnt = TRUE)[[1L]]
      rebuilt <- GCEstimjmv:::.gcePrecisionState(state$precisionBasis, which, c(1, 1, 1))
      testthat::expect_equal(rebuilt$data,
        GCEstimjmv:::.gceDiagnosticState(original, which)$data)
    }

    # A fresh analysis with only persisted cache mimics a reopened analysis
    # whose images/tables must be regenerated. Any refit is a test failure.
    testthat::local_mocked_bindings(
      cv.lmgce = function(...) stop("Unexpected CV refit"),
      ridgetrace = function(...) stop("Unexpected ridge refit"),
      .package = "GCEstim")
    unavailable <- cache_analysis(kind, data, list(trueCoef = "garbage", plot5 = TRUE, plot7 = TRUE))
    option <- unavailable$options$options$supportSignalVectorN
    option$value <- 1L
    option <- unavailable$options$options$twostepsN
    option$value <- 0L
    unavailable$results$fitCache$setState(state)
    unavailable$run()
    testthat::expect_false(".fitWarning" %in% names(unavailable$results))
    testthat::expect_false(unavailable$results$plot5$visible)
    testthat::expect_false(unavailable$results$plot7$visible)

    # Disabled input, invalid active input, and corrected input all reuse the fit.
    for (input in c("garbage", "1, 2", "1, Inf, 1", "1, 1,", "")) {
      for (selected in c(FALSE, TRUE)) {
        invalid <- cache_analysis(kind, data, list(trueCoef = input,
          plot1 = TRUE, plot2 = TRUE, plot5 = selected, plot7 = selected))
        invalid$results$fitCache$setState(state)
        invalid$run()
        testthat::expect_equal(invalid$results$coefficients$asDF, first$results$coefficients$asDF)
        for (key in c("plot1", "plot2"))
          testthat::expect_gt(nrow(GCEstimjmv:::.gcePlotImage(invalid$results[[key]], state)$state$data), 0L)
        testthat::expect_null(invalid$results$plot5$state)
        testthat::expect_null(invalid$results$plot7$state)
        testthat::expect_identical(".fitWarning" %in% names(invalid$results), selected)
        if (selected)
          testthat::expect_match(invalid$results[[".fitWarning"]]$content, "Precision-error plots skipped")
      }
    }
    changedTruth <- cache_analysis(kind, data, list(trueCoef = "2, 2, 2", plot5 = TRUE, plot7 = TRUE))
    changedTruth$results$fitCache$setState(state)
    changedTruth$run()
    testthat::expect_equal(changedTruth$results$coefficients$asDF, first$results$coefficients$asDF)
    testthat::expect_equal(changedTruth$results$plot5$state$data,
      GCEstimjmv:::.gcePrecisionState(state$precisionBasis, 5L, c(2, 2, 2))$data)
    table_names <- c("modelSummary", "signalSupports", "coefficients",
                     "signalProbabilities", "cvResults")
    retained <- lapply(table_names, function(name) first$results[[name]]$asDF)
    # Loaded tables can survive option changes and .init(). Reusing them must
    # preserve row counts and values rather than appending duplicate rows.
    first$.__enclos_env__$private$.init()
    first$run()
    for (i in seq_along(table_names))
      testthat::expect_equal(first$results[[table_names[i]]]$asDF, retained[[i]])
    for (enabled in c(FALSE, TRUE)) {
      replay <- cache_analysis(kind, data, list(
        plot1 = TRUE, plot2 = TRUE, plot3 = TRUE, plot4 = TRUE,
        plot5 = TRUE, plot6 = TRUE, plot7 = TRUE, plotCV = TRUE,
        OLS = enabled, NormEnt = enabled))
      for (key in c("saveFitted", "saveResiduals", "saveW")) {
        option <- replay$options$options[[key]]
        option$value <- list(value = TRUE)
      }
      replay$results$fitCache$.__enclos_env__$super$fromProtoBuf(state_message(state),
                                          character(), character())
      replay$run()
      testthat::expect_equal(replay$results$coefficients$asDF,
                             first$results$coefficients$asDF)
      testthat::expect_equal(replay$results$signalProbabilities$asDF,
                             first$results$signalProbabilities$asDF)
      testthat::expect_equal(replay$results$cvResults$asDF,
                             first$results$cvResults$asDF)
      testthat::expect_identical(replay$results$fitCache$state, state)
      for (key in names(state$plots)) {
        testthat::expect_identical(replay$results[[key]]$state$key, key)
        testthat::expect_null(replay$results[[key]]$state$data)
      }
      # Rendering still works after jamovi serializes both cache and image state.
      device <- tempfile(fileext = ".pdf")
      grDevices::pdf(device)
      tryCatch({
        for (key in c(names(state$plots), "plot5", "plot7")) {
          renderer <- replay$.__enclos_env__$private[[paste0(".", key)]]
          testthat::expect_true(renderer(replay$results[[key]], ggplot2::theme_minimal(),
                                       list(color = c("blue", "green"))))
        }
      }, finally = { grDevices::dev.off(); unlink(device) })
      testthat::expect_identical(GCEstimjmv:::.gcePlotImage(replay$results$plot2, state)$state$showNEP, enabled)
      testthat::expect_identical(replay$results$plot5$state$showNEP, enabled)
      testthat::expect_identical(
        any(is.finite(GCEstimjmv:::.gcePlotImage(replay$results$plot1, state)$state$data$ols)), enabled)
      testthat::expect_identical(!is.null(GCEstimjmv:::.gcePlotImage(replay$results$plot3, state)$state$ols), enabled)
      for (key in c("plot2", "plot5", "plot6", "plot7"))
        testthat::expect_identical(any(vapply(
          GCEstimjmv:::.gcePlotImage(replay$results[[key]], state)$state$overlays,
          function(x) "yintercept" %in% names(x), logical(1))), enabled)
      for (key in c("saveFitted", "saveResiduals", "saveW")) {
        output <- replay$results[[key]]$.__enclos_env__$private
        expected <- switch(key, saveFitted = cached$fit$best$fitted.values,
          saveResiduals = cached$fit$best$residuals, saveW = cached$fit$best$w)
        testthat::expect_equal(unlist(output$.rowNums),
          which(stats::complete.cases(data[, c("y", "X001", "X002")])))
        testthat::expect_equal(unname(do.call(cbind, output$.values)),
                               unname(as.matrix(expected)))
      }
    }
  })
})

testthat::test_that("diagnostic state preserves entropy uncertainty", {
  data <- data.frame(support = 1:2, error = c(2, 3), nep = c(.4, .6),
                     nep.LL = c(.3, .5), nep.UL = c(.5, .7))
  plot <- ggplot2::ggplot(data)
  state <- GCEstimjmv:::.gceDiagnosticState(plot, 2L)
  testthat::expect_equal(state$data, data)
})
