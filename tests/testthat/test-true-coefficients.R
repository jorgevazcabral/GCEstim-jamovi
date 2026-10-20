# Stop at the estimator boundary: these tests exercise option validation,
# independently of whether a candidate support converges in GCEstim.
for (kind in c("TARWGCE", "TASWGCE")) local({
  kind <- kind
  testthat::test_that(paste(kind, "ignores disabled stale coefficients after model changes"), {
    data <- utils::read.csv(testthat::test_path("_fixtures", "dataThesis.csv"))
    ns <- asNamespace("GCEstimjmv")
    testthat::local_mocked_bindings(cv.lmgce = function(...) stop("estimator reached"),
                                   .package = "GCEstim")
    for (change in list(list(factors = "X005"), list(intercept = FALSE),
      list(supportSignalVectorN = 1, plot5 = TRUE, plot7 = TRUE))) {
      for (input in c("1, 1, 1", "not numeric")) {
        args <- modifyList(list(dep = "y", covs = c("X001", "X002"), run = TRUE,
          M = "3", J = "3", weight = "0.5", cvNfolds = 2,
          supportSignalVectorN = 2, twostepsN = 0, bootstrap = FALSE,
          plot1 = TRUE, plot2 = FALSE, plot5 = FALSE, plot7 = FALSE,
          plotCV = FALSE, trueCoef = input, seed = 230676), change)
        if (kind == "TARWGCE") args$ridgeLambdaN <- 100
        options <- do.call(get(paste0(kind, "Options"), ns)$new, args)
        analysis <- get(paste0(kind, "Class"), ns)$new(options = options, data = data)
        testthat::expect_error(analysis$run(), "estimator reached")
        testthat::expect_false(analysis$results$plot5$visible)
        testthat::expect_false(analysis$results$plot7$visible)
      }
    }
  })
})

testthat::test_that("true coefficient validation handles factor and intercept dimensions", {
  data <- utils::read.csv(testthat::test_path("_fixtures", "dataThesis.csv"))
  validate <- GCEstimjmv:::.gceTrueCoefficients
  testthat::expect_length(validate("1, 1, 1, 1, 1", y ~ X001 + X002 + X005, data)$values, 5L)
  testthat::expect_match(validate("1, 1, 1", y ~ X001 + X002 + X005, data)$warning, "5 values")
  testthat::expect_match(validate("1, 1, 1", y ~ X001 + X002 - 1, data)$warning, "2 values")
  testthat::expect_match(validate("1", y ~ missing, data)$warning, "skipped")
})
