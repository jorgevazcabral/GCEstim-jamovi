wgce_cached_analysis <- function(data, changes = list()) {
  ns <- asNamespace("GCEstimjmv")
  options <- do.call(get("WGCEOptions", ns)$new, modifyList(list(
    dep = "y", covs = c("X001", "X002"), factors = NULL,
    run = TRUE, plot1 = FALSE, bootstrap = FALSE, seed = 230676,
    bootB = 10), changes))
  get("WGCEClass", ns)$new(options = options, data = data)
}

wgce_plain_state <- function(x) {
  if (inherits(x, c("lmgce", "ggplot")) || is.environment(x) || is.function(x)) return(FALSE)
  if (is.list(x) && !all(vapply(x, wgce_plain_state, logical(1)))) return(FALSE)
  attrs <- attributes(x)
  is.null(attrs) || all(vapply(attrs, wgce_plain_state, logical(1)))
}

for (method in c("z", "percentile", "basic")) local({
  method <- method
  testthat::test_that(paste("W-GCE replays", method, "intervals, plots and saved columns"), {
    data <- utils::read.csv(testthat::test_path("_fixtures", "dataThesis.csv"))
    data$X001[c(2, 9)] <- NA_real_
    changes <- list(bootstrap = method != "z")
    if (method != "z") changes$bootCIMethod <- method
    original <- GCEstim::lmgce
    testthat::local_mocked_bindings(lmgce = function(...) {
      warning("fit warning retained")
      original(...)
    }, .package = "GCEstim")
    first <- wgce_cached_analysis(data, changes)
    first$run()
    state <- unserialize(serialize(first$results$fitCache$state, NULL))
    testthat::expect_identical(state$version, 4L)
    testthat::expect_true(wgce_plain_state(state))
    testthat::expect_true(all(is.finite(state$ci)))
    testthat::expect_true(any(is.finite(state$plots$plot1$data$ols)))
    testthat::expect_match(first$results[[".fitWarning"]]$content, "fit warning retained")
    testthat::local_mocked_bindings(lmgce = function(...) stop("Unexpected fit"), .package = "GCEstim")
    testthat::local_mocked_bindings(confint = function(...) stop("Unexpected CI"), .package = "stats")
    tables <- c("modelSummary", "signalSupports", "coefficients", "signalProbabilities")
    retained <- lapply(tables, function(key) first$results[[key]]$asDF)
    first$results$signalProbabilities$getColumn("p1")$setTitle("p_1")
    first$.__enclos_env__$private$.init()
    first$run()
    testthat::expect_identical(first$results$signalProbabilities$getColumn("p1")$title, "p (point 1)")
    for (i in seq_along(tables))
      testthat::expect_equal(first$results[[tables[i]]]$asDF, retained[[i]])
    for (ols in c(FALSE, TRUE, FALSE)) {
      replay <- wgce_cached_analysis(data, modifyList(changes, list(plot1 = TRUE, OLS = ols)))
      for (key in c("saveFitted", "saveResiduals", "saveW")) {
        option <- replay$options$options[[key]]
        option$value <- list(value = TRUE)
      }
      # Match reopening: restore only the compressed persisted cache state.
      replay$results$fitCache$.__enclos_env__$super$fromProtoBuf(state_message(state), character(), character())
      replay$run()
      testthat::expect_identical(replay$results$fitCache$state, state)
      for (i in seq_along(tables))
        testthat::expect_equal(replay$results[[tables[i]]]$asDF, retained[[i]])
      testthat::expect_match(replay$results[[".fitWarning"]]$content, "fit warning retained")
      testthat::expect_match(replay$results[[".fitTime"]]$content, "Previous fit reused")
      image <- GCEstimjmv:::.gcePlotImage(replay$results$plot1, state)
      testthat::expect_equal(image$state$data$lower, state$plots$plot1$data$lower)
      testthat::expect_equal(image$state$data$upper, state$plots$plot1$data$upper)
      testthat::expect_identical(any(is.finite(image$state$data$ols)), ols)
      testthat::expect_null(replay$results$plot1$state$data)
      device <- tempfile(fileext = ".pdf")
      grDevices::pdf(device)
      tryCatch(testthat::expect_true(replay$.__enclos_env__$private$.plot1(
        replay$results$plot1, ggplot2::theme_minimal(), list(color = "blue"))),
        finally = { grDevices::dev.off(); unlink(device) })
      for (key in c("saveFitted", "saveResiduals", "saveW")) {
        output <- replay$results[[key]]$.__enclos_env__$private
        expected <- switch(key, saveFitted = state$fit$fitted.values,
          saveResiduals = state$fit$residuals, saveW = state$fit$w)
        testthat::expect_equal(unlist(output$.rowNums),
          which(stats::complete.cases(data[, c("y", "X001", "X002")])))
        testthat::expect_equal(unname(do.call(cbind, output$.values)), unname(as.matrix(expected)))
      }
    }
  })
})

testthat::test_that("W-GCE cache depends on every fit and CI option, never display or save choices", {
  data <- utils::read.csv(testthat::test_path("_fixtures", "dataThesis.csv"))
  expected <- c("data", "run", "dep", "covs", "factors", "intercept",
    "supportSignalLowers", "supportSignalUppers", "M", "J", "weight",
    "uniformSignalPrior", "signalPriorWeights", "uniformNoisePrior", "noisePriorWeights",
    "supportNoise3sig", "supportNoiseMin", "supportNoiseMax", "method",
    "bootConfLevel", "bootstrap", "bootB", "bootMethod", "bootCIMethod", "seed")
  first <- wgce_cached_analysis(data)
  testthat::expect_setequal(first$results$fitCache$.__enclos_env__$private$.clearWith, expected)
  saved <- state_message(list(sentinel = TRUE))
  for (change in c(expected, "plot1", "OLS", "saveFitted", "saveResiduals", "saveW")) {
    replay <- wgce_cached_analysis(data)
    replay$results$fitCache$.__enclos_env__$super$fromProtoBuf(saved, change, character())
    testthat::expect_identical(is.null(replay$results$fitCache$state), change %in% expected, info = change)
  }
})

testthat::test_that("W-GCE missing or incompatible caches refit once", {
  data <- utils::read.csv(testthat::test_path("_fixtures", "dataThesis.csv"))
  fit <- GCEstim::lmgce
  ci <- stats::confint
  fits <- intervals <- 0L
  testthat::local_mocked_bindings(lmgce = function(...) { fits <<- fits + 1L; fit(...) }, .package = "GCEstim")
  testthat::local_mocked_bindings(confint = function(...) { intervals <<- intervals + 1L; ci(...) }, .package = "stats")
  for (old in list(NULL, list(version = 3L))) {
    analysis <- wgce_cached_analysis(data)
    analysis$results$fitCache$setState(old)
    analysis$run()
    analysis$run()
  }
  testthat::expect_identical(fits, 2L)
  testthat::expect_identical(intervals, 2L)
})


testthat::test_that("fit, CI and prior changes refit a restored W-GCE cache", {
  data <- utils::read.csv(testthat::test_path("_fixtures", "dataThesis.csv"))
  first <- wgce_cached_analysis(data)
  first$run()
  state <- first$results$fitCache$state
  fit <- GCEstim::lmgce
  ci <- stats::confint
  fits <- intervals <- 0L
  testthat::local_mocked_bindings(lmgce = function(...) { fits <<- fits + 1L; fit(...) }, .package = "GCEstim")
  testthat::local_mocked_bindings(confint = function(...) { intervals <<- intervals + 1L; ci(...) }, .package = "stats")
  changes <- list(list(seed = 123), list(bootConfLevel = .9),
    list(uniformSignalPrior = FALSE, signalPriorWeights = "0.1, 0.2, 0.7"),
    list(uniformNoisePrior = FALSE, noisePriorWeights = "0.1, 0.2, 0.7"))
  for (change in changes) {
    analysis <- wgce_cached_analysis(data, change)
    analysis$results$fitCache$.__enclos_env__$super$fromProtoBuf(state_message(state), names(change), character())
    testthat::expect_null(analysis$results$fitCache$state)
    analysis$run()
    analysis$run()
  }
  testthat::expect_identical(fits, 4L)
  testthat::expect_identical(intervals, 4L)
})

testthat::test_that("adaptive support point hints use the same prefix and retain examples", {
  for (kind in c("tarwgce", "taswgce")) {
    lines <- readLines(testthat::test_path("..", "..", "jamovi", paste0(kind, ".u.yaml")))
    hints <- grep("odd integers at least 3", lines, value = TRUE, ignore.case = TRUE)
    testthat::expect_length(hints, 2L)
    testthat::expect_true(all(grepl("One or more odd integers at least 3; ", hints, fixed = TRUE)))
    testthat::expect_true(any(grepl("e.g., 3, 9, 15.", hints, fixed = TRUE)))
    testthat::expect_true(any(grepl("e.g., 3, 5, 7.", hints, fixed = TRUE)))
  }
})
