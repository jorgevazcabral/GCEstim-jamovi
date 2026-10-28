low_data <- function() read.csv(testthat::test_path("_fixtures", "dataThesis.csv"))
low_analysis <- function(changes = list()) {
  ns <- asNamespace("GCEstimjmv")
  options <- do.call(get("RidgeTraceOptions", ns)$new, modifyList(list(dep = "y", covs = c("X001", "X002"),
    run = TRUE, lambdaN = 8L, cvNfolds = 3L, plotSample = "show", plotRanges = TRUE,
    plotSensitivity = TRUE, plotFolds = TRUE, plotStability = TRUE), changes))
  get("RidgeTraceClass", ns)$new(options = options, data = low_data())
}
low_keys <- c("plotPaths", "plotCV", "plotSample", "plotRanges", "plotSensitivity", "plotFolds", "plotStability")

testthat::test_that("LOW dependencies invalidate exactly the plots that read each option", {
  expected <- list(plotPaths = c("selection", "manualIndex", "axisScale", "showIntercept", "extrema", "facet", "displayTerms"),
    plotCV = c("axisScale", "seBars"), plotSample = c("selection", "manualIndex", "axisScale"),
    plotRanges = c("selection", "manualIndex", "showIntercept", "displayTerms"),
    plotSensitivity = c("selection", "manualIndex", "axisScale", "showIntercept", "displayTerms"),
    plotFolds = c("selection", "manualIndex", "axisScale", "foldMean"),
    plotStability = c("selection", "manualIndex", "axisScale", "showIntercept", "displayTerms"))
  estimation <- c("data", "run", "dep", "covs", "factors", "intercept", "standardize", "penalizeIntercept", "customGrid",
    "lambdaValues", "lambdaMin", "lambdaMax", "lambdaN", "cv", "cvNfolds", "seed", "errorMeasure")
  presentation <- unique(unlist(expected))
  for (key in low_keys) {
    image <- low_analysis()$results[[key]]
    testthat::expect_setequal(image$.__enclos_env__$private$.clearWith, c(estimation, expected[[key]], key))
    for (option in c(estimation, presentation, low_keys)) {
      replay <- low_analysis()$results[[key]]
      replay$.__enclos_env__$super$fromProtoBuf(state_message(list(previous = TRUE)), option, character())
      testthat::expect_identical(is.null(replay$state), option %in% c(estimation, expected[[key]], key), info = paste(key, option))
    }
  }
})

testthat::test_that("LOW unrelated options retain actual plot states and avoid rebuilding them", {
  a <- low_analysis(); a$run(); cache <- a$results$fitCache$state
  original <- GCEstimjmv:::.gceRTPlot; built <- character()
  testthat::local_mocked_bindings(.gceRTPlot = function(cache, state, ...) {
    built <<- c(built, state$key); original(cache, state, ...)
  }, .package = "GCEstimjmv")
  testthat::local_mocked_bindings(ridgetrace = function(...) stop("Unnecessary refit"), .package = "GCEstim")
  cases <- list(seBars = list(value = FALSE, keys = "plotCV"), foldMean = list(value = FALSE, keys = "plotFolds"),
    extrema = list(value = TRUE, keys = "plotPaths"), facet = list(value = TRUE, keys = "plotPaths"),
    displayTerms = list(value = "X001", keys = c("plotPaths", "plotRanges", "plotSensitivity", "plotStability")),
    showIntercept = list(value = FALSE, keys = c("plotPaths", "plotRanges", "plotSensitivity", "plotStability")),
    axisScale = list(value = "linear", keys = setdiff(low_keys, "plotRanges")))
  for (key in names(cases)) {
    # Reset to a known presentation before checking one change.
    b <- low_analysis(); b$results$fitCache$setState(cache); b$run(); built <- character()
    before <- lapply(low_keys, function(k) b$results[[k]]$state); names(before) <- low_keys
    option <- b$options$option(key); option$value <- cases[[key]]$value
    b$run()
    testthat::expect_setequal(built, cases[[key]]$keys)
    for (other in setdiff(low_keys, cases[[key]]$keys)) testthat::expect_identical(b$results[[other]]$state, before[[other]], info = paste(key, other))
    testthat::expect_identical(b$results$fitCache$state, cache)
  }
  built <- character(); a$run(); testthat::expect_length(built, 0L)
  beforeCV <- a$results$plotCV$state
  option <- a$options$option("selection"); option$value <- "manual"
  option <- a$options$option("manualIndex"); option$value <- 3L
  a$run()
  testthat::expect_setequal(built, setdiff(low_keys, "plotCV"))
  testthat::expect_identical(a$results$plotCV$state, beforeCV)
  for (key in setdiff(low_keys, "plotCV")) testthat::expect_identical(a$results[[key]]$state$selection$index, 3L)
})

testthat::test_that("LOW effective standardization preserves preferences through intercept transitions", {
  engine <- GCEstim::ridgetrace; calls <- list()
  testthat::local_mocked_bindings(ridgetrace = function(...) {
    args <- list(...); calls[[length(calls)+1L]] <<- args; do.call(engine, args)
  }, .package = "GCEstim")
  a <- low_analysis()
  for (intercept in c(TRUE, FALSE, TRUE, FALSE, TRUE)) {
    option <- a$options$option("intercept"); option$value <- intercept; a$run()
    testthat::expect_true(a$options$standardize)
    testthat::expect_identical(tail(calls, 1)[[1]]$standardize, intercept)
    testthat::expect_identical(a$results$modelSummary$asDF$value[a$results$modelSummary$asDF$measure == "Standardization"], if (intercept) "Yes" else "Not applicable")
    direct <- engine(if (intercept) y ~ X001 + X002 else y ~ X001 + X002 - 1, low_data(),
      standardize = intercept, lambda.n = 8, cv.nfolds = 3)
    testthat::expect_identical(a$results$fitCache$state$fit$coefficients, unname(direct$coef.lambda) |> `rownames<-`(rownames(direct$coef.lambda)))
  }
  testthat::expect_length(calls, 5L)
  b <- low_analysis(list(intercept = FALSE)); b$run(); before <- b$results$fitCache$state; count <- length(calls)
  option <- b$options$option("standardize"); option$value <- FALSE; b$run()
  testthat::expect_identical(b$results$fitCache$state, before)
  testthat::expect_length(calls, count)
  ui <- yaml::read_yaml(file.path(testthat::test_path(), "../../jamovi/ridgetrace.u.yaml"))
  testthat::expect_identical(Filter(function(x) identical(x$name, "standardize"), ui$children)[[1]]$enable, "(intercept)")
})

testthat::test_that("LOW CV fallback preserves stored 1-SE and uses effective rule everywhere", {
  a <- low_analysis(list(selection = "1se"))
  for (cv in c(TRUE, FALSE, TRUE, FALSE, TRUE)) {
    option <- a$options$option("cv"); option$value <- cv; a$run()
    f <- a$results$fitCache$state$fit
    expected <- GCEstimjmv:::.gceRTSelection(f, if (cv) "1se" else "min", 1, cv)
    testthat::expect_identical(a$options$selection, "1se")
    testthat::expect_identical(a$results$plotPaths$state$selection$rule, if (cv) "1se" else "min")
    testthat::expect_identical(a$results$plotPaths$state$selection$index, expected$index)
    testthat::expect_identical(a$results$coefficients$asDF$estimate, as.numeric(f$coefficients[,expected$index]))
    testthat::expect_match(a$results$selectedLambda$asDF$rule, if (cv) "1-SE" else "descriptive")
    testthat::expect_identical(a$results$selectedLambda$asDF$selected, f$lambda[expected$index])
    if (!cv) {
      testthat::expect_match(a$results[[".method"]]$content, "CV is disabled; minimum in-sample error")
      testthat::expect_identical(a$results$plotSample$state$selection$rule, "min")
    }
  }
})

testthat::test_that("LOW numeric lambda results retain precision and concise coefficient titles", {
  a <- low_analysis(list(customGrid = TRUE, lambdaValues = "0.123456789012345, 2.987654321098765", selection = "manual", manualIndex = 2L, showLambdaTable = TRUE))
  a$run(); f <- a$results$fitCache$state$fit; table <- a$results$selectedLambda$asDF
  for (key in c("minimum", "maximum", "selected")) {
    testthat::expect_type(table[[key]], "double")
    testthat::expect_equal(a$results$selectedLambda$getColumn(key)$type, "number")
  }
  testthat::expect_identical(table$selected, f$lambda[2])
  testthat::expect_identical(table$minimum, min(f$lambda))
  testthat::expect_identical(table$maximum, max(f$lambda))
  testthat::expect_identical(table$rule, "Manual (grid position)")
  testthat::expect_match(a$results$coefficients$title, "2.99", fixed = TRUE)
  testthat::expect_identical(a$results$lambdaResults$asDF$lambda, f$lambda)
  testthat::expect_false(any(c("Minimum lambda", "Maximum lambda", "Selected lambda") %in% a$results$modelSummary$asDF$measure))
  testthat::expect_identical(a$results$lambdaResults$asDF$selected, c(NA_character_, "Yes"))
})

testthat::test_that("LOW captions are absent and caveats remain in table notes", {
  a <- low_analysis(); a$run()
  for (key in low_keys) {
    p <- GCEstimjmv:::.gceRTPlot(a$results$fitCache$state, a$results[[key]]$state, NULL, list())
    testthat::expect_null(p$labels$caption, info = key)
    testthat::expect_silent(ggplot2::ggplot_build(p))
  }
  testthat::expect_match(a$results[[".method"]]$content, "complete dataset")
  testthat::expect_false(grepl("confidence intervals", a$results[[".method"]]$content))
  # Source-level footnote checks complement the rendered table APIs.
  source <- paste(readLines(file.path(testthat::test_path(), "../../R/ridgetrace.b.R")), collapse = "\n")
  testthat::expect_match(source, "Ranges describe the evaluated lambda grid; they are not confidence intervals.", fixed = TRUE)
  testthat::expect_false(grepl("selectedLambda$setNote", source, fixed = TRUE))
  testthat::expect_false(grepl("coefficients$setNote", source, fixed = TRUE))
  testthat::expect_false(grepl("grid-dependent finite-difference sensitivity", source, fixed = TRUE))
  testthat::expect_match(source, "CV summaries require all folds to be finite", fixed = TRUE)
})

testthat::test_that("LOW result and option titles follow the approved labels", {
  a <- low_analysis()
  titles <- c("Ridge Coefficient Paths", "Cross-Validation Error", "In-Sample Error", "Coefficient Ranges",
    "Local Coefficient Sensitivity", "Fold-Specific CV Errors", "Coefficient Stability Across CV Folds")
  testthat::expect_identical(vapply(low_keys, function(k) a$results[[k]]$title, "", USE.NAMES = FALSE), titles)
  definitions <- yaml::yaml.load(paste(readLines(file.path(testthat::test_path(), "../../jamovi/ridgetrace.a.yaml"), encoding = "UTF-8"), collapse = "\n"))
  options <- setNames(definitions$options, vapply(definitions$options, `[[`, "", "name"))
  testthat::expect_identical(options$seBars$title, "CV mean \u00B1 1 SE")
  testthat::expect_identical(options$selection$options[[3]]$title, "Manual (grid position)")
  for (key in setdiff(low_keys, "plotSample")) testthat::expect_identical(options[[key]]$title, a$results[[key]]$title)
  files <- c("R/ridgetrace-state.R", "R/ridgetrace.b.R")
  source <- paste(unlist(lapply(file.path(testthat::test_path(), "../..", files), readLines)), collapse = "\n")
  testthat::expect_false(grepl(".removeOption", source, fixed = TRUE))
  testthat::expect_false(grepl(".gceRTClean", source, fixed = TRUE))
})

testthat::test_that("LOW genuine v1.1.0 numerical and image states reopen without refitting", {
  cache <- readRDS(testthat::test_path("_fixtures", "ridge-1.1.0-fitCache.rds"))
  image <- readRDS(testthat::test_path("_fixtures", "ridge-1.1.0-plotPaths.rds"))
  testthat::expect_identical(cache$version, 1L)
  testthat::expect_null(cache$settings)
  testthat::expect_false(any(c("plotHeatmap", "plotSelected", "normalize") %in% names(image$options)))
  data <- low_data(); names(data)[2] <- "X002 (f)"; data$X001[40] <- NA_real_; data$X003[14] <- NA_real_
  ns <- asNamespace("GCEstimjmv")
  options <- get("RidgeTraceOptions", ns)$new(dep = "y", covs = c("X001", "X002 (f)", "X003", "X004"),
    run = TRUE, displayTerms = "X002 (f)")
  a <- get("RidgeTraceClass", ns)$new(options = options, data = data)
  a$results$fitCache$.__enclos_env__$super$fromProtoBuf(state_message(cache), character(), character())
  a$results$plotPaths$.__enclos_env__$super$fromProtoBuf(state_message(image), character(), character())
  testthat::local_mocked_bindings(ridgetrace = function(...) stop("Released cache refitted"), .package = "GCEstim")
  a$run()
  testthat::expect_identical(a$results$fitCache$state, cache)
  testthat::expect_identical(a$results$plotPaths$state$terms, "`X002 (f)`")
  testthat::expect_identical(a$results$plotPaths$state$presentationVersion, 2L)
  testthat::expect_identical(a$results$selectedLambda$asDF$selected, cache$fit$lambda[59L])
  testthat::expect_identical(a$results$coefficients$asDF$estimate, as.numeric(cache$fit$coefficients[,59L]))
  testthat::expect_match(a$results[[".fitTime"]]$content, "Previous fit reused")
})

testthat::test_that("LOW entropy recovery failures use WARNING notices while retaining valid rows", {
  ns <- asNamespace("GCEstimjmv")
  rows <- data.frame(support = c(1,2), predictor = "X", estimate = c(1,2), nepk = c(.2,.3))
  state <- list(which = 4L, coefficientTerms = "X", data = rows[2,,drop=FALSE], overlays = list())
  for (kind in c("TARWGCE", "TASWGCE")) {
    options <- get(paste0(kind, "Options"), ns)$new(plot4 = TRUE)
    results <- get(paste0(kind, "Results"), ns)$new(options = options)
    GCEstimjmv:::.gceEntropyWarning(results, options, list(plot4 = state))
    testthat::expect_true(results[[".entropyWarning"]]$visible)
    testthat::expect_match(results[[".entropyWarning"]]$content, "Only available observations")
    recovered <- GCEstimjmv:::.gceEntropyState(state)
    testthat::expect_identical(recovered$data, state$data)
    GCEstimjmv:::.gceEntropyWarning(results, options, list(plot4 = state, plot3 = list(which = 3L, data = rows)))
    testthat::expect_false(results[[".entropyWarning"]]$visible)
    testthat::expect_identical(GCEstimjmv:::.gceEntropyState(state, list(which = 3L, data = rows))$data, rows)
  }
})
