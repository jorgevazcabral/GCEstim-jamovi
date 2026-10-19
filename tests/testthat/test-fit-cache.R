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
    first <- cache_analysis(kind, data)
    first$run()
    state <- unserialize(serialize(first$results$fitCache$state, NULL))
    cached <- GCEstimjmv:::.gceUnpackFitCache(state$payload)
    testthat::expect_identical(cached$version, 3L)
    testthat::expect_true(cache_plain(cached))
    testthat::expect_true(all(c("plot1", paste0("plot", 2:7), "plotCV") %in%
                                names(cached$plots)))
    testthat::expect_true(all(vapply(cached$plots, function(x)
      is.data.frame(x$data) && nrow(x$data) > 0, logical(1))))
    clearWith <- first$results$fitCache$.__enclos_env__$private$.clearWith
    for (change in c("OLS", "NormEnt", "plot1", "plot5", "saveFitted",
                     "saveResiduals", "saveW", "cvTableRows")) {
      testthat::expect_false(change %in% clearWith)
    }
    for (change in c("data", "M", "seed", "bootConfLevel", "trueCoef")) {
      testthat::expect_true(change %in% clearWith)
    }

    # A fresh analysis with only persisted cache mimics a reopened analysis
    # whose images/tables must be regenerated. Any refit is a test failure.
    testthat::local_mocked_bindings(
      cv.lmgce = function(...) stop("Unexpected CV refit"),
      ridgetrace = function(...) stop("Unexpected ridge refit"),
      .package = "GCEstim")
    for (enabled in c(FALSE, TRUE)) {
      replay <- cache_analysis(kind, data, list(
        plot1 = TRUE, plot2 = TRUE, plot3 = TRUE, plot4 = TRUE,
        plot5 = TRUE, plot6 = TRUE, plot7 = TRUE, plotCV = TRUE,
        OLS = enabled, NormEnt = enabled))
      for (key in c("saveFitted", "saveResiduals", "saveW")) {
        option <- replay$options$options[[key]]
        option$value <- list(value = TRUE)
      }
      replay$results$fitCache$setState(state)
      replay$run()
      testthat::expect_equal(replay$results$coefficients$asDF,
                             first$results$coefficients$asDF)
      testthat::expect_equal(replay$results$signalProbabilities$asDF,
                             first$results$signalProbabilities$asDF)
      testthat::expect_equal(replay$results$cvResults$asDF,
                             first$results$cvResults$asDF)
      testthat::expect_identical(replay$results$fitCache$state, state)
      testthat::expect_identical(replay$results$plot2$state$showNEP, enabled)
      testthat::expect_identical(replay$results$plot5$state$showNEP, enabled)
      testthat::expect_identical(
        any(is.finite(replay$results$plot1$state$data$ols)), enabled)
      testthat::expect_identical(!is.null(replay$results$plot3$state$ols), enabled)
      for (key in c("plot2", "plot5", "plot6", "plot7"))
        testthat::expect_identical(any(vapply(
          replay$results[[key]]$state$overlays,
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
