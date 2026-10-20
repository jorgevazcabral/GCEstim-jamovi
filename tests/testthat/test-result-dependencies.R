testthat::test_that("TARW results retain precise invalidation dependencies", {
  ns <- asNamespace("GCEstimjmv")
  tarw <- get("TARWGCEResults", ns)$new(
    options = get("TARWGCEOptions", ns)$new())
  tasw <- get("TASWGCEResults", ns)$new(
    options = get("TASWGCEOptions", ns)$new())
  clear_with <- function(result) result$.__enclos_env__$private$.clearWith
  ridge_options <- c("ridgeSignalSupport", "noiseSupportMethod",
                     "ridgeLambdaMin", "ridgeLambdaMax", "ridgeLambdaN")
  for (name in c("coefficients", "signalProbabilities")) {
    testthat::expect_setequal(clear_with(tarw[[name]]),
      c(clear_with(tasw[[name]]), ridge_options))
  }
  for (name in c("cvResults", "plotCV")) {
    testthat::expect_false("OLS" %in% clear_with(tarw[[name]]))
    testthat::expect_setequal(clear_with(tarw[[name]]),
      c(clear_with(tasw[[name]]), ridge_options))
  }
  testthat::expect_setequal(clear_with(tarw$plotRidge), c(
    "data", "run", "dep", "covs", "factors", "intercept",
    "ridgeLambdaMin", "ridgeLambdaMax", "ridgeLambdaN", "errorMeasure",
    "seed", "plotRidge"))
})

testthat::test_that("W-GCE invalidates affected results for all four prior options", {
  ns <- asNamespace("GCEstimjmv")
  results <- get("WGCEResults", ns)$new(options = get("WGCEOptions", ns)$new())
  priors <- c("uniformSignalPrior", "signalPriorWeights",
              "uniformNoisePrior", "noisePriorWeights")
  for (name in c("modelSummary", "signalSupports", "coefficients",
                 "plot1", "signalProbabilities")) {
    dependencies <- results[[name]]$.__enclos_env__$private$.clearWith
    testthat::expect_true(all(priors %in% dependencies), info = name)
    results[[name]]$setState(list(previous = TRUE))
    saved <- state_message(results[[name]]$state)
    for (prior in priors) {
      fresh <- get("WGCEResults", ns)$new(options = get("WGCEOptions", ns)$new())
      fresh[[name]]$.__enclos_env__$super$fromProtoBuf(saved, prior, character())
      testthat::expect_null(fresh[[name]]$state, info = paste(name, prior))
    }
    if (name != "plot1")
      testthat::expect_false(any(c("OLS", "plot1") %in% dependencies))
  }
})

testthat::test_that("changing W-GCE priors changes estimates and rejects invalid weights", {
  data <- utils::read.csv(testthat::test_path("_fixtures", "dataThesis.csv"))
  fit <- function(...) GCEstimjmv::WGCE(data = data, dep = y,
    covs = vars(X001, X002), factors = NULL, run = TRUE, plot1 = FALSE, bootstrap = FALSE, ...)
  baseline <- fit()
  signal <- fit(uniformSignalPrior = FALSE, signalPriorWeights = "0.1, 0.2, 0.7")
  noise <- fit(uniformNoisePrior = FALSE, noisePriorWeights = "0.1, 0.2, 0.7")
  for (changed in list(signal, noise)) {
    testthat::expect_false(isTRUE(all.equal(changed$coefficients$asDF$estimate,
                                          baseline$coefficients$asDF$estimate)))
  }
  testthat::expect_error(fit(uniformSignalPrior = FALSE, signalPriorWeights = "bad"))
  testthat::expect_error(fit(uniformNoisePrior = FALSE, noisePriorWeights = "bad"))
})
