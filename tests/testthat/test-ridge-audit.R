audit_analysis <- function(data, changes = list()) {
  ns <- asNamespace("GCEstimjmv")
  options <- do.call(get("RidgeTraceOptions", ns)$new, modifyList(list(
    dep = "y", covs = c("X001", "X002"), run = TRUE, lambdaN = 8L,
    cvNfolds = 3L, plotSample = "show", plotRanges = TRUE,
    plotSensitivity = TRUE, plotFolds = TRUE, plotStability = TRUE), changes))
  get("RidgeTraceClass", ns)$new(options = options, data = data)
}
audit_data <- function() read.csv(testthat::test_path("_fixtures", "dataThesis.csv"))

testthat::test_that("audit colour maps cover 5, 7, 10 and 15 actual model coefficients", {
  actual <- jmvcore:::getGlobalTheme("default", "jmv")
  set.seed(128)
  for (n in c(5L, 7L, 10L, 15L)) {
    d <- as.data.frame(matrix(rnorm(80 * (n-1L)), nrow = 80))
    d$y <- rowSums(d) + rnorm(80)
    a <- audit_analysis(d, list(covs = setdiff(names(d), "y"))); a$run()
    cache <- a$results$fitCache$state; f <- cache$fit; terms <- rownames(f$coefficients)
    colours <- GCEstimjmv:::.gceRTColours(f, actual$theme)
    testthat::expect_length(colours, n)
    testthat::expect_equal(length(unique(colours)), n)
    # CIE Lab distance records a reproducible minimum separation check.
    lab <- grDevices::convertColor(t(grDevices::col2rgb(colours))/255, "sRGB", "Lab")
    testthat::expect_gt(min(stats::dist(lab)), 5)
    reordered <- f; reordered$coefficients <- f$coefficients[rev(terms), , drop = FALSE]
    testthat::expect_identical(GCEstimjmv:::.gceRTColours(reordered, actual$theme), colours)
    testthat::expect_equal(length(unique(GCEstimjmv:::.gceRTColours(f, list(color = c("blue", "green"))))), n)
    for (key in c("plotPaths", "plotRanges", "plotSensitivity", "plotStability")) {
      for (hide in c(FALSE, TRUE)) for (filter in c("", terms[2])) {
        state <- a$results[[key]]$state; state$terms <- NULL
        state$options$showIntercept <- !hide; state$options$displayTerms <- filter
        p <- suppressMessages(GCEstimjmv:::.gceRTPlot(cache, state, actual$ggtheme, actual$theme))
        built <- ggplot2::ggplot_build(p); visible <- unique(as.character(p$data$term))
        testthat::expect_identical(built$plot$scales$get_scales("colour")$map(visible), unname(colours[visible]))
        testthat::expect_identical("(Intercept)" %in% visible, !hide && filter == "")
      }
    }
  }
})

testthat::test_that("audit filters resolve displayed, internal, factor and quoted comma names", {
  terms <- c("(Intercept)", "`income per person`", "`group name`level B", "`region, east`", "regular")
  labels <- GCEstimjmv:::.gceCoefficientLabels(terms)
  testthat::expect_identical(GCEstimjmv:::.gceRTLabels(terms), labels)
  for (i in seq_along(terms)) for (name in c(terms[i], labels[i])) {
    text <- paste0('"', name, '"')
    testthat::expect_identical(GCEstimjmv:::.gceRTResolveTerms(terms, text), terms[i])
  }
  testthat::expect_identical(GCEstimjmv:::.gceRTResolveTerms(terms, '"region, east", income per person'), terms[c(2,4)])
  for (text in c("unknown", '"unclosed')) {
    error <- tryCatch(GCEstimjmv:::.gceRTResolveTerms(terms, text), error = identity)
    testthat::expect_s3_class(error, "error")
    for (label in labels) testthat::expect_match(conditionMessage(error), label, fixed = TRUE)
  }
  collisions <- c("`same name`", "same name", "same name [coefficient] 1")
  labels <- GCEstimjmv:::.gceCoefficientLabels(collisions)
  testthat::expect_false(anyDuplicated(labels) > 0)
  testthat::expect_error(GCEstimjmv:::.gceRTResolveTerms(collisions, "same name"), "Ambiguous.*Valid displayed")
  testthat::expect_identical(GCEstimjmv:::.gceRTResolveTerms(collisions, '`same name`'), collisions[1])
  testthat::expect_identical(GCEstimjmv:::.gceRTResolveTerms(collisions, labels[2]), collisions[2])
  d <- audit_data(); names(d)[2] <- "income, person"; names(d)[3] <- "group name"
  d[[3]] <- factor(d[[3]] > median(d[[3]]), labels = c("level A", "level B"))
  a <- audit_analysis(d, list(covs = "income, person", factors = "group name")); a$run()
  cache <- a$results$fitCache$state
  for (label in a$results$coefficients$asDF$term) {
    b <- audit_analysis(d, list(covs = "income, person", factors = "group name", displayTerms = paste0('"', label, '"')))
    b$results$fitCache$setState(cache); b$run()
    testthat::expect_length(b$results$plotPaths$state$terms, 1L)
    testthat::expect_identical(b$results$fitCache$state, cache)
    testthat::expect_identical(b$results$coefficients$asDF, a$results$coefficients$asDF)
  }
})

testthat::test_that("audit invalid filters reject once before fitting and recover using saved cache", {
  a <- audit_analysis(audit_data()); a$run(); cache <- a$results$fitCache$state
  testthat::local_mocked_bindings(ridgetrace = function(...) stop("Unnecessary refit"), .package = "GCEstim")
  invalid <- audit_analysis(audit_data(), list(displayTerms = "unknown"))
  invalid$results$fitCache$setState(cache)
  testthat::expect_warning(testthat::expect_error(invalid$run(), "Unknown coefficient names.*Valid displayed"), NA)
  testthat::expect_true(invalid$results$plotPaths$visible)
  testthat::expect_identical(invalid$results$fitCache$state, cache)
  testthat::expect_false(".fitWarning" %in% names(invalid$results))
  fresh <- audit_analysis(audit_data(), list(displayTerms = "unknown"))
  testthat::expect_error(fresh$run(), "Unknown coefficient names")
  option <- invalid$options$option("displayTerms")
  option$value <- "X001"
  invalid$run()
  testthat::expect_identical(invalid$results$plotPaths$state$terms, "X001")
  testthat::expect_identical(invalid$results$fitCache$state, cache)
  testthat::expect_identical(invalid$results$coefficients$asDF, a$results$coefficients$asDF)
  testthat::expect_identical(invalid$results$plotPaths$state$selection, a$results$plotPaths$state$selection)
})

testthat::test_that("audit unexpected build failures preserve visible plot errors and other results", {
  a <- audit_analysis(audit_data()); a$run(); cache <- a$results$fitCache$state
  original <- GCEstimjmv:::.gceRTPlot
  testthat::local_mocked_bindings(.gceRTPlot = function(cache, state, ...) {
    if (state$key == "plotPaths") stop("injected build failure")
    original(cache, state, ...)
  }, .package = "GCEstimjmv")
  b <- audit_analysis(audit_data()); b$results$fitCache$setState(cache)
  testthat::expect_warning(b$run(), NA)
  image <- b$results$plotPaths
  testthat::expect_true(image$visible)
  testthat::expect_identical(image$status, "error")
  testthat::expect_match(image$state$error, paste0(image$title, ": injected build failure"), fixed = TRUE)
  testthat::expect_false(grepl("plotPaths", image$state$error, fixed = TRUE))
  testthat::expect_error(GCEstimjmv:::.gceRTRender(image, cache, NULL, list()), image$title, fixed = TRUE)
  testthat::expect_identical(b$results$coefficients$asDF, a$results$coefficients$asDF)
  testthat::expect_identical(b$results$fitCache$state, cache)
  for (key in c("plotCV", "plotSample", "plotRanges", "plotSensitivity", "plotFolds", "plotStability")) {
    testthat::expect_true(b$results[[key]]$visible)
    testthat::expect_identical(b$results[[key]]$status, "complete")
  }
  testthat::expect_false(".fitWarning" %in% names(b$results))
})

testthat::test_that("audit unexpected render errors propagate with the user-visible title", {
  a <- audit_analysis(audit_data()); a$run(); cache <- a$results$fitCache$state
  testthat::local_mocked_bindings(.gceRTPlot = function(...) stop("injected render failure"), .package = "GCEstimjmv")
  image <- a$results$plotRanges
  testthat::expect_error(GCEstimjmv:::.gceRTRender(image, cache, NULL, list()), paste0(image$title, ": injected render failure"), fixed = TRUE)
  testthat::expect_true(image$visible)
  testthat::expect_identical(a$results$fitCache$state, cache)
  testthat::expect_true(a$results$plotCV$visible)
  testthat::expect_identical(image$status, "error")
})
