# Capture renderer output without opening a graphics device or changing namespaces.
graphics_plot <- function(name, state, global = jmvcore:::getGlobalTheme("default", "jmv")) {
  fun <- get(name, asNamespace("GCEstimjmv"))
  capture <- new.env(parent = environment(fun))
  capture$print <- function(x) { capture$plot <- x; invisible(x) }
  environment(fun) <- capture
  testthat::expect_true(suppressMessages(fun(list(state = state), global$ggtheme, global$theme)))
  list(plot = capture$plot, built = ggplot2::ggplot_build(capture$plot))
}
graphics_analysis <- function(kind, data, changes = list()) {
  args <- list(dep = "y", covs = c("Z signal", "A signal"), factors = "F group", run = TRUE, plot1 = TRUE)
  if (kind != "WGCE") args <- c(args, list(M = "3", J = "3", weight = "0.5", cvNfolds = 5L,
    supportSignalVectorN = 5L, twostepsN = 0L, plot2 = TRUE, plot3 = TRUE, plot4 = TRUE, plotCV = FALSE))
  if (kind == "TARWGCE") args <- c(args, list(ridgeLambdaN = 20L, plotRidge = TRUE))
  ns <- asNamespace("GCEstimjmv")
  options <- do.call(get(paste0(kind, "Options"), ns)$new, modifyList(args, changes))
  get(paste0(kind, "Class"), ns)$new(options = options, data = data)
}
graphics_data <- function() {
  d <- read.csv(testthat::test_path("_fixtures", "dataThesis.csv"))
  names(d)[names(d) == "X001"] <- "Z signal"
  names(d)[names(d) == "X002"] <- "A signal"
  names(d)[names(d) == "X005"] <- "F group"
  d[["F group"]] <- factor(d[["F group"]]); d
}
graphics_map <- function(p, ids) unname(p$built$plot$scales$get_scales("colour")$map(ids))
graphics_points <- function(p) Filter(function(x) inherits(x$geom, "GeomPoint") &&
  is.data.frame(x$data) && "predictor" %in% names(x$data), p$plot$layers)

testthat::test_that("existing-analysis labels are presentation only and collision safe", {
  ids <- c("(Intercept)", "`A space`", "A space", "A space [coefficient] 1", "`F group`level B", "a! + (b)", "a`b")
  labels <- GCEstimjmv:::.gceCoefficientLabels(ids)
  testthat::expect_identical(names(labels), ids)
  testthat::expect_false(anyDuplicated(labels) > 0L)
  testthat::expect_identical(unname(labels[c(2, 5, 6, 7)]), c("A space", "F grouplevel B", "a! + (b)", "a`b"))
  testthat::expect_identical(labels, GCEstimjmv:::.gceRTLabels(ids))
})

testthat::test_that("named coefficient scales survive actual themes, data ordering and subsets", {
  ids <- c("(Intercept)", "`Z signal`", "`A signal`", "FB", "FC", "D", "E")
  s <- list(coefficientTerms = ids, data = data.frame(term = ids, estimate = 1:7,
    lower = 0:6, upper = 2:8, ols = 2:8), caption = "Estimate")
  for (palette in c("jmv", "hadley", "spss")) {
    g <- jmvcore:::getGlobalTheme("default", palette)
    expected <- jmvcore::colorPalette(7L, palette, "color")
    for (rows in list(seq_along(ids), c(1,3,2,7,4,6,5), 2:7, c(2,5))) {
      state <- s; state$data <- state$data[rows, ]
      p <- graphics_plot(".gceRenderEstimates", unserialize(serialize(state, NULL)), g)
      testthat::expect_identical(graphics_map(p, ids[rows]), expected[rows])
      testthat::expect_identical(as.character(p$plot$data$term), ids[rows])
      testthat::expect_identical(p$built$data[[3]]$colour, expected[rows])
      testthat::expect_identical(p$built$data[[2]]$colour, expected[rows])
      testthat::expect_identical(p$built$data[[4]]$colour, expected[rows])
      testthat::expect_identical(as.vector(p$built$layout$panel_params[[1]]$y$get_labels(), "character"),
        unname(GCEstimjmv:::.gceCoefficientLabels(ids)[ids[rows]]))
    }
  }
})

for (kind in c("WGCE", "TASWGCE", "TARWGCE")) local({
  kind <- kind
  testthat::test_that(paste(kind, "preserves values and named identities through legacy and current caches"), {
    d <- graphics_data(); ns <- asNamespace("GCEstimjmv"); full <- NULL
    if (kind == "WGCE") {
      original <- GCEstim::lmgce
      testthat::local_mocked_bindings(lmgce = function(...) { full <<- original(...); full }, .package = "GCEstim")
    } else {
      original <- GCEstim::cv.lmgce
      testthat::local_mocked_bindings(cv.lmgce = function(...) { full <<- original(...); full }, .package = "GCEstim")
    }
    a <- graphics_analysis(kind, d); a$run(); cache <- a$results$fitCache$state
    best <- if (kind == "WGCE") full else full$best
    ids <- names(stats::coef(best)); labels <- GCEstimjmv:::.gceCoefficientLabels(ids)
    expected <- jmvcore::colorPalette(max(6L, length(ids)), "jmv", "color")[seq_along(ids)]
    testthat::expect_identical(rownames(cache$coefficients), ids)
    testthat::expect_identical(cache$version, 4L)
    testthat::expect_equal(a$results$coefficients$asDF$estimate, unname(stats::coef(best)))
    ci <- stats::confint(best)
    testthat::expect_equal(a$results$coefficients$asDF$ciLower, unname(ci[ids, 1]))
    testthat::expect_equal(a$results$coefficients$asDF$ciUpper, unname(ci[ids, 2]))
    for (key in c("coefficients", "signalSupports", "signalProbabilities"))
      testthat::expect_identical(a$results[[key]]$asDF$term, unname(labels))
    testthat::expect_identical(cache$plots$plot1$data$term, ids)
    estimate <- graphics_plot(".gceRenderEstimates", GCEstimjmv:::.gcePlotImage(a$results$plot1, cache)$state)
    testthat::expect_identical(graphics_map(estimate, ids), expected)
    testthat::expect_equal(estimate$built$data[[3]]$x, a$results$coefficients$asDF$estimate)
    testthat::expect_equal(estimate$built$data[[2]]$xmin, a$results$coefficients$asDF$ciLower)
    testthat::expect_equal(estimate$built$data[[2]]$xmax, a$results$coefficients$asDF$ciUpper)
    zero <- estimate$built$data[[1]]
    testthat::expect_true(all(zero$xintercept == 0 & zero$linetype == "dashed" & zero$colour == "grey50"))
    if (kind != "WGCE") for (key in c("plot3", "plot4")) {
      state <- GCEstimjmv:::.gcePlotImage(a$results[[key]], cache)$state
      p <- graphics_plot(".gceRenderDiagnostic", state)
      testthat::expect_identical(graphics_map(p, ids), expected)
      selected <- graphics_points(p)[[1]]$data
      testthat::expect_identical(selected$predictor, ids)
      testthat::expect_equal(selected$y, if (key == "plot3") unname(stats::coef(best)) else
        unname(GCEstim::NormEnt(best, model = FALSE)))
      testthat::expect_identical(graphics_map(p, selected$predictor), expected)
      index <- which(vapply(p$plot$layers, function(x) inherits(x$geom, "GeomPoint"), logical(1)))
      testthat::expect_identical(p$built$data[[index]]$colour, expected)
      testthat::expect_identical(unname(p$built$plot$scales$get_scales("colour")$get_labels(ids)), unname(labels))
    }
    ridge <- if (kind == "TARWGCE") a$results$plotRidge$state else NULL
    if (!is.null(ridge)) {
      p <- graphics_plot(".gceRenderRidge", ridge)
      testthat::expect_identical(graphics_map(p, ids), expected)
      testthat::expect_identical(p$plot$theme$legend.position, "bottom")
      testthat::expect_s3_class(p$plot$theme$legend.title, "element_blank")
      testthat::expect_false(any(vapply(p$plot$layers, function(x) inherits(x$geom, "GeomVline"), logical(1))))
      points <- Filter(function(x) inherits(x$geom, "GeomPoint"), p$plot$layers)
      testthat::expect_identical(unname(vapply(points, function(x) x$aes_params$shape, numeric(1))), c(24,25))
      testthat::expect_identical(unname(vapply(points, function(x) x$aes_params$size, numeric(1))), c(1,1))
      for (id in ids) {
        path <- ridge$data[ridge$data$variable == id, ]; extrema <- ridge$extrema[ridge$extrema$variable == id, ]
        testthat::expect_equal(extrema$max_y, max(path$coefficient))
        testthat::expect_equal(extrema$min_y, min(path$coefficient))
        testthat::expect_equal(extrema$max_x, path$lambda[which.max(path$coefficient)])
        testthat::expect_equal(extrema$min_x, path$lambda[which.min(path$coefficient)])
      }
      testthat::expect_identical(p$built$data[[2]]$colour, expected)
      testthat::expect_identical(p$built$data[[3]]$colour, expected)
    }
    testthat::local_mocked_bindings(lmgce = function(...) stop("Unexpected estimation"),
      cv.lmgce = function(...) stop("Unexpected CV"), ridgetrace = function(...) stop("Unexpected ridge estimation"), .package = "GCEstim")
    # Also refresh already-populated old table labels, without changing numerical cells.
    before <- a$results$coefficients$asDF; keys <- unlist(a$results$coefficients$rowKeys)
    for (i in seq_along(ids)) a$results$coefficients$setRow(rowKey = keys[i], values = list(term = ids[i]))
    a$run()
    testthat::expect_equal(a$results$coefficients$asDF, before)
    testthat::expect_identical(a$results$fitCache$state, cache)
    for (legacy in c(FALSE, TRUE)) {
      saved <- unserialize(serialize(cache, NULL)); savedRidge <- ridge
      if (legacy) {
        for (key in names(saved$plots)) {
          saved$plots[[key]]$coefficientTerms <- NULL
          saved$plots[[key]]$overlays <- lapply(saved$plots[[key]]$overlays, function(x) { x$predictor <- NULL; x })
        }
        if (!is.null(savedRidge)) savedRidge$coefficientTerms <- NULL
      }
      replay <- graphics_analysis(kind, d, list(OLS = TRUE))
      replay$results$fitCache$.__enclos_env__$super$fromProtoBuf(state_message(saved), character(), character())
      if (!is.null(savedRidge)) replay$results$plotRidge$.__enclos_env__$super$fromProtoBuf(state_message(savedRidge), character(), character())
      replay$run()
      testthat::expect_identical(replay$results$fitCache$state, saved)
      testthat::expect_equal(replay$results$coefficients$asDF, before)
      testthat::expect_match(replay$results[[".fitTime"]]$content, "Previous fit reused")
      for (key in if (kind == "WGCE") "plot1" else c("plot1", "plot3", "plot4")) {
        state <- GCEstimjmv:::.gcePlotImage(replay$results[[key]], saved)$state
        p <- graphics_plot(if (key == "plot1") ".gceRenderEstimates" else ".gceRenderDiagnostic", state)
        testthat::expect_identical(graphics_map(p, ids), expected)
        if (key != "plot1") {
          points <- graphics_points(p)[[1]]$data
          testthat::expect_identical(points$predictor, ids)
          originalPoint <- Filter(function(x) all(c("x", "y") %in% names(x)), cache$plots[[key]]$overlays)[[1]]
          testthat::expect_identical(points$x, originalPoint$x)
          testthat::expect_identical(points$y, originalPoint$y)
        }
      }
      if (!is.null(savedRidge)) {
        p <- graphics_plot(".gceRenderRidge", replay$results$plotRidge$state)
        testthat::expect_identical(graphics_map(p, ids), expected)
        testthat::expect_identical(p$plot$theme$legend.position, "bottom")
      }
    }
  })
})

testthat::test_that("support references and non-coefficient palettes are unchanged", {
  s <- list(which = 3L, data = data.frame(support = rep(1:3, each = 2), predictor = rep(c("Z", "A"),3), estimate = 1:6),
    overlays = list(), ols = NULL, supportLines = list(min = 1, oneSE = 2, manual = 3))
  p <- graphics_plot(".gceRenderDiagnostic", s)
  refs <- Filter(function(x) "xintercept" %in% names(x), p$built$data)
  testthat::expect_identical(vapply(refs, function(x) x$linetype[1], character(1)), c("dotted", "dashed", "dotted"))
  testthat::expect_identical(vapply(refs, function(x) x$colour[1], character(1)), c("black", "black", "green4"))
  testthat::expect_identical(vapply(refs, function(x) x$xintercept[1], numeric(1)), c(1,2,3))
  s <- list(which = 2L, data = data.frame(support = 1:3, error = c(1,2,4), nep = c(.2,.4,.6)),
    overlays = list(), supportLines = list(), showNEP = TRUE)
  p <- graphics_plot(".gceRenderDiagnostic", s)
  g <- jmvcore:::getGlobalTheme("default", "jmv")
  testthat::expect_true(all(p$built$data[[1]]$colour == g$theme$color[1]))
  testthat::expect_true(all(p$built$data[[2]]$colour == g$theme$color[2]))
  s <- list(data = data.frame(signalPoints = rep(c(3,5),4), noisePoints = rep(c(3,5), each=2, times=2),
    weight = rep(c(.3,.5),each=4), cvMean = 1:8), measure = "RMSE")
  p <- graphics_plot(".gceRenderCV", s)
  for (panel in unique(p$built$data[[1]]$PANEL))
    testthat::expect_identical(unique(p$built$data[[1]]$colour[p$built$data[[1]]$PANEL == panel]), jmvcore::colorPalette(2L,"jmv","color"))
})
