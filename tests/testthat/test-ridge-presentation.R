testthat::test_that("ridge presentation maps preserve internal identities and disambiguate labels", {
  terms <- c("(Intercept)", "X001", "`X002 (f)`", "X003", "X004", "X002 (f)", "X002 (f) [coefficient] 1")
  labels <- GCEstimjmv:::.gceRTLabels(terms)
  testthat::expect_identical(names(labels), terms)
  testthat::expect_identical(unname(labels[3]), "X002 (f)")
  testthat::expect_false(anyDuplicated(labels) > 0L)
  testthat::expect_identical(unname(GCEstimjmv:::.gceRTLabels("`factor name`level B")), "factor namelevel B")
})

testthat::test_that("ridge colours, labels, exact values and line types survive filtering and reopening", {
  ns <- asNamespace("GCEstimjmv")
  make <- function(changes = list()) {
    o <- do.call(get("RidgeTraceOptions", ns)$new, modifyList(list(dep = "y", covs = c("X001", "X002"),
      run = TRUE, lambdaN = 8L, cvNfolds = 3L, plotSample = "show", plotRanges = TRUE,
      plotSensitivity = TRUE, plotFolds = TRUE, plotStability = TRUE), changes))
    get("RidgeTraceClass", ns)$new(options = o, data = read.csv(testthat::test_path("_fixtures", "dataThesis.csv")))
  }
  a <- make(); a$run(); cache <- a$results$fitCache$state
  # Use internal names with a quoted factor label without altering numerical estimates.
  terms <- rownames(cache$fit$coefficients); terms[3] <- "`X002 (f)`"
  rownames(cache$fit$coefficients) <- terms
  cache$fit$foldCoefficients <- lapply(cache$fit$foldCoefficients, function(b) { rownames(b) <- terms; b })
  theme <- list(color = c("lightblue", "grey", "orange", "green", "red"))
  expected <- GCEstimjmv:::.gceRTColours(cache$fit, theme)
  testthat::expect_identical(expected, GCEstimjmv:::.gceRTColours(unserialize(serialize(cache$fit, NULL)), theme))
  testthat::local_mocked_bindings(ridgetrace = function(...) stop("Unexpected refit"), .package = "GCEstim")
  for (rule in c("min", "1se", "manual")) for (hide in c(FALSE, TRUE)) {
    b <- make(list(selection = rule, manualIndex = 3L, showIntercept = !hide))
    b$results$fitCache$.__enclos_env__$super$fromProtoBuf(state_message(cache), character(), character()); b$run()
    testthat::expect_identical(b$results$fitCache$state, cache)
    testthat::expect_identical(b$results$coefficients$asDF$term, unname(GCEstimjmv:::.gceRTLabels(terms)))
    for (key in c("plotPaths", "plotRanges", "plotSensitivity", "plotStability", "plotSample", "plotFolds")) {
      state <- unserialize(serialize(b$results[[key]]$state, NULL))
      for (filter in c("", terms[2])) {
        state$options$displayTerms <- filter
        state$terms <- NULL # exercise legacy-image reconstruction with changed options
        p <- suppressMessages(GCEstimjmv:::.gceRTPlot(cache, state, list(ggplot2::theme_minimal(), ggplot2::scale_colour_manual(values = rev(theme$color))), theme))
        built <- ggplot2::ggplot_build(p)
        if (key %in% c("plotPaths", "plotRanges", "plotSensitivity", "plotStability")) {
          visible <- unique(as.character(p$data$term))
          testthat::expect_identical(built$plot$scales$get_scales("colour")$map(visible), unname(expected[visible]))
          testthat::expect_identical(unname(built$plot$scales$get_scales("colour")$get_labels(terms)), unname(GCEstimjmv:::.gceRTLabels(terms)))
        }
        if (key == "plotRanges") {
          positions <- match(as.character(p$data$term), terms)
          testthat::expect_equal(p$data$estimate, unname(cache$fit$coefficients[positions, state$selection$index]))
          testthat::expect_equal(p$data$minimum, cache$fit$minimum[positions])
          testthat::expect_equal(p$data$maximum, cache$fit$maximum[positions])
          testthat::expect_equal(built$data[[3]]$x, p$data$estimate)
        } else {
          last <- tail(built$data, 1)[[1]]
          testthat::expect_true(all(last$linetype == switch(rule, min = "dotted", `1se` = "dashed", manual = "solid")))
          testthat::expect_true(all(last$xintercept == log10(cache$fit$lambda[state$selection$index])))
        }
      }
    }
  }
  # CV-disabled automatic and manual selections retain the same conventions.
  for (rule in c("min", "manual")) {
    state <- a$results$plotPaths$state; state$selection$rule <- rule; state$terms <- NULL
    noCV <- cache; noCV$fit$cvErrors <- NULL; noCV$fit$foldCoefficients <- NULL
    p <- GCEstimjmv:::.gceRTPlot(noCV, state, NULL, theme)
    testthat::expect_identical(tail(p$layers, 1)[[1]]$aes_params$linetype, if (rule == "min") "dotted" else "solid")
  }
})

testthat::test_that("reported example retains its exact grid position and selected estimates", {
  d <- read.csv(testthat::test_path("_fixtures", "dataThesis.csv"))
  d$X001[40] <- NA_real_; d$X003[14] <- NA_real_
  names(d)[names(d) == "X002"] <- "X002 (f)"
  a <- GCEstimjmv::RidgeTrace(data=d, dep="y", covs=c("X001","X002 (f)","X003","X004"), factors=NULL,
    run=TRUE, plotRanges=TRUE, facet=TRUE)
  f <- a$fitCache$state$fit; s <- a$plotPaths$state$selection
  testthat::expect_identical(a$plotPaths$state$presentationVersion, 2L)
  testthat::expect_identical(s$index, 59L)
  testthat::expect_equal(f$lambda[s$index], 3.2745491628777286, tolerance=1e-14)
  expected <- c(-4.045055104986524, .6404993249794027, -13.999103646322222, -12.150484293055118, -4.9108547996093055)
  testthat::expect_equal(a$coefficients$asDF$estimate, expected, tolerance=1e-12)
  testthat::expect_equal(a$coefficientRanges$asDF$selected, expected, tolerance=1e-12)
  testthat::expect_identical(a$coefficients$asDF$term, c("(Intercept)","X001","X002 (f)","X003","X004"))
  theme <- list(color=c("lightblue","grey","orange","green","red"))
  colours <- GCEstimjmv:::.gceRTColours(f, theme)
  testthat::expect_identical(unname(colours[rownames(f$coefficients)]), theme$color)
  p <- GCEstimjmv:::.gceRTPlot(a$fitCache$state, a$plotPaths$state, NULL, theme)
  built <- ggplot2::ggplot_build(p)
  testthat::expect_identical(built$plot$scales$get_scales("colour")$map(rownames(f$coefficients)), theme$color)
  actualTheme <- jmvcore:::getGlobalTheme("default", "jmv")
  p <- suppressMessages(GCEstimjmv:::.gceRTPlot(a$fitCache$state, a$plotRanges$state, actualTheme$ggtheme, actualTheme$theme))
  built <- ggplot2::ggplot_build(p)
  testthat::expect_identical(built$data[[3]]$colour, jmvcore::colorPalette(5L, "jmv", "color"))
  testthat::expect_identical(rownames(f$coefficients)[3], "`X002 (f)`")
})
