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
