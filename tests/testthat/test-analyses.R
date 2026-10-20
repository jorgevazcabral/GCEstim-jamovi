thesis_data <- function() {
  # utils::data() reads package .csv datasets with a semicolon separator;
  # this copy of the bundled dataThesis.csv is comma-delimited.
  data <- utils::read.csv(
    testthat::test_path("_fixtures", "dataThesis.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  testthat::expect_s3_class(data, "data.frame")
  testthat::expect_gt(nrow(data), 2L)
  testthat::expect_true(all(c("y", "X001", "X002") %in% names(data)))
  testthat::expect_true(is.numeric(data$y))
  data
}

expect_filled_table <- function(results, name, columns) {
  table <- results[[name]]$asDF
  testthat::expect_s3_class(table, "data.frame")
  testthat::expect_gt(nrow(table), 0L)
  testthat::expect_true(all(columns %in% names(table)))
  table
}

expect_successful_fit <- function(results, adaptive = FALSE) {
  # Warnings are inserted dynamically as a Notice named .fitWarning.
  testthat::expect_false(".fitWarning" %in% names(results))

  summary <- expect_filled_table(results, "modelSummary", c("measure", "value"))
  testthat::expect_true(any(nzchar(as.character(summary$value))))
  testthat::expect_identical(as.character(summary$value[summary$measure == "Confidence interval type"]),
                             "Asymptotic normal")
  supports <- expect_filled_table(results, "signalSupports", c("term", "lower", "upper"))
  testthat::expect_true(any(is.finite(supports$lower) & is.finite(supports$upper)))

  coefficients <- expect_filled_table(
    results, "coefficients", c("term", "estimate", "ciLower", "ciUpper")
  )
  testthat::expect_true(any(is.finite(coefficients$estimate)))
  # In particular, catch the default (non-bootstrap) hasCI regression.
  testthat::expect_true(any(is.finite(coefficients$ciLower) &
                            is.finite(coefficients$ciUpper)))

  probabilities <- expect_filled_table(results, "signalProbabilities", "term")
  p_columns <- grep("^p[0-9]+$", names(probabilities), value = TRUE)
  testthat::expect_gt(length(p_columns), 0L)
  for (j in seq_along(p_columns))
    testthat::expect_identical(results$signalProbabilities$getColumn(p_columns[j])$title,
                               paste0("p (point ", j, ")"))
  testthat::expect_true(any(is.finite(as.numeric(probabilities[[p_columns[1L]]]))))

  if (adaptive) {
    cv <- expect_filled_table(
      results, "cvResults", c("signalPoints", "noisePoints", "cvMean")
    )
    testthat::expect_true(any(is.finite(cv$cvMean)))
  }

  invisible(coefficients)
}

testthat::test_that("W-GCE fills its tables and default confidence intervals", {
  data <- thesis_data()
  results <- GCEstimjmv::WGCE(
    data = data, dep = y, covs = vars(X001, X002), factors = NULL,
    run = TRUE, bootstrap = FALSE, seed = 230676
  )
  expect_successful_fit(results)
})

testthat::test_that("W-GCE accepts a predictor name containing a space", {
  data <- thesis_data()
  names(data)[names(data) == "X001"] <- "Signal X001"
  results <- GCEstimjmv::WGCE(
    data = data, dep = y, covs = c("Signal X001", "X002"), factors = NULL,
    run = TRUE, bootstrap = FALSE, seed = 230676
  )
  coefficients <- expect_successful_fit(results)
  # The model summary may retain R's backticks around a non-syntactic name.
  term_labels <- gsub("`", "", as.character(coefficients$term), fixed = TRUE)
  testthat::expect_true("Signal X001" %in% term_labels)
})

testthat::test_that("TARW-GCE fills its model and CV tables", {
  results <- GCEstimjmv::TARWGCE(
    data = thesis_data(), dep = y, covs = vars(X001, X002), factors = NULL,
    run = TRUE, M = "3, 5", J = "3, 5", weight = "0.3, 0.5", cvNfolds = 5,
    cvTableRows = 8, supportSignalVectorN = 5, ridgeLambdaN = 100,
    twostepsN = 0, bootstrap = FALSE, plot1 = FALSE, plot2 = FALSE,
    plotCV = FALSE, plotRidge = FALSE, seed = 230676
  )
  expect_successful_fit(results, adaptive = TRUE)
})

testthat::test_that("TASW-GCE fills its model and CV tables", {
  results <- GCEstimjmv::TASWGCE(
    data = thesis_data(), dep = y, covs = vars(X001, X002), factors = NULL,
    run = TRUE, M = "3", J = "3", weight = "0.5", cvNfolds = 2,
    cvTableRows = 1, supportSignalVectorN = 2, twostepsN = 0,
    bootstrap = FALSE, plot1 = FALSE, plot2 = FALSE, plotCV = FALSE,
    seed = 230676
  )
  expect_successful_fit(results, adaptive = TRUE)
})
