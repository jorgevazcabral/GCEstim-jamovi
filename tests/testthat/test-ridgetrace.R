rt_plain <- function(x) {
  if (is.environment(x) || is.function(x) || inherits(x, c("ridgetrace", "ggplot"))) return(FALSE)
  if (is.list(x) && !all(vapply(x, rt_plain, logical(1)))) return(FALSE)
  is.null(attributes(x)) || all(vapply(attributes(x), rt_plain, logical(1)))
}
rt_analysis <- function(data, changes = list()) {
  ns <- asNamespace("GCEstimjmv")
  options <- do.call(get("RidgeTraceOptions", ns)$new, modifyList(list(
    dep = "y", covs = c("X001", "X002"), factors = NULL, run = TRUE,
    lambdaN = 8L, cvNfolds = 3L, showLambdaTable = TRUE), changes))
  get("RidgeTraceClass", ns)$new(options = options, data = data)
}
rt_data <- function() utils::read.csv(testthat::test_path("_fixtures", "dataThesis.csv"))

testthat::test_that("ridge tables and orientation agree with GCEstim across model settings", {
  d <- rt_data(); d$X001[2] <- NA_real_
  for (settings in list(list(), list(cv = FALSE), list(standardize = FALSE),
      list(standardize = FALSE, penalizeIntercept = FALSE),
      list(standardize = FALSE, intercept = FALSE), list(factors = "X005"))) {
    a <- rt_analysis(d, settings); a$run(); o <- a$options
    complete <- d[complete.cases(d[, c("y", "X001", "X002", unlist(o$factors))]), ]
    form <- if (length(o$factors)) y ~ X001 + X002 + X005 else if (!o$intercept) y ~ X001 + X002 - 1 else y ~ X001 + X002
    direct <- GCEstim::ridgetrace(form, complete, lambda.n = 8, cv.nfolds = 3,
      cv = o$cv, standardize = o$standardize, penalize.intercept = o$penalizeIntercept)
    f <- a$results$fitCache$state$fit
    testthat::expect_equal(f$lambda, direct$lambda)
    testthat::expect_equal(unname(f$coefficients), unname(direct$coef.lambda))
    testthat::expect_identical(rownames(f$coefficients), rownames(direct$coef.lambda))
    testthat::expect_equal(f$sample, direct$error.lambda)
    s <- GCEstimjmv:::.gceRTSelection(f, "min", .001, o$cv)
    testthat::expect_equal(a$results$coefficients$asDF$estimate, as.numeric(direct$coef.lambda[, s$index]))
    testthat::expect_equal(a$results$coefficientRanges$asDF$minimum, as.numeric(direct$min.coef))
    testthat::expect_equal(a$results$coefficientRanges$asDF$maximum, as.numeric(direct$max.coef))
    testthat::expect_equal(a$results$coefficientRanges$asDF$maxAbs, as.numeric(direct$max.abs.coef))
    if (o$cv) {
      testthat::expect_equal(s$mean, unname(colMeans(direct$error.lambda.cv)))
      testthat::expect_equal(s$se, unname(apply(direct$error.lambda.cv, 2, sd) / sqrt(3)))
      testthat::expect_match(a$results[[".method"]]$content, if (o$standardize) "complete dataset" else "Fold-based SE")
    } else {
      testthat::expect_null(f$cvErrors)
      testthat::expect_true(a$results$plotSample$visible)
      testthat::expect_false(a$results$plotCV$visible)
    }
  }
})

testthat::test_that("ridge custom grid uses numeric positions despite rounded name collisions", {
  a <- rt_analysis(rt_data(), list(customGrid = TRUE, lambdaValues = "1, 0.0000000012, 0.0000000011"))
  a$run(); f <- a$results$fitCache$state$fit
  testthat::expect_identical(f$lambda, c(1.1e-9, 1.2e-9, 1))
  testthat::expect_equal(ncol(f$coefficients), 3L)
  direct <- GCEstim::ridgetrace(y ~ X001 + X002, rt_data(), lambda = f$lambda, cv.nfolds = 3)
  testthat::expect_equal(unname(f$coefficients), unname(direct$coef.lambda))
  for (text in c("", "1", "1,1", "0,1", "-1,1", "Inf,1", "1,", "bad,1"))
    testthat::expect_error(GCEstimjmv:::.gceRTLambda(TRUE, text, .001, 1000, 100), "Custom lambda")
})

testthat::test_that("ridge minimum, 1-SE and manual selection have defined semantics", {
  f <- list(lambda = c(.1, 1, 10, 100), sample = c(1,2,3,4),
    cvErrors = rbind(c(2,1,1.1,3), c(2,1.4,1.3,3), c(2,.6,1.2,3)))
  minimum <- GCEstimjmv:::.gceRTSelection(f, "min", 1, TRUE)
  one <- GCEstimjmv:::.gceRTSelection(f, "1se", 1, TRUE)
  testthat::expect_equal(minimum$index, 2L)
  testthat::expect_equal(one$index, 3L)
  testthat::expect_equal(one$se[2], sd(c(1,1.4,.6))/sqrt(3))
  testthat::expect_equal(GCEstimjmv:::.gceRTSelection(f, "manual", 100, TRUE)$index, 4L)
  testthat::expect_error(GCEstimjmv:::.gceRTSelection(f, "manual", 2, TRUE), "evaluated")
  testthat::expect_error(GCEstimjmv:::.gceRTSelection(f, "1se", 1, FALSE), "requires cross-validation")
  close <- f; close$lambda <- c(1, 1 + .Machine$double.eps, 2, 3)
  testthat::expect_equal(GCEstimjmv:::.gceRTSelection(close, "manual", close$lambda[2], TRUE)$index, 2L)
  f$cvErrors[1,2] <- Inf
  testthat::expect_true(is.na(GCEstimjmv:::.gceRTSelection(f, "min", 1, TRUE)$mean[2]))
  f$cvErrors[,] <- NaN
  testthat::expect_true(is.na(GCEstimjmv:::.gceRTSelection(f, "min", 1, TRUE)$index))
  testthat::expect_equal(GCEstimjmv:::.gceRTSelection(f, "manual", 1, TRUE)$index, 2L)
})

testthat::test_that("ridge unsupported settings and unexpected failures are informative", {
  d <- rt_data()
  testthat::expect_error({ a <- rt_analysis(d, list(intercept = FALSE)); a$run() }, "implementation limitation")
  testthat::expect_error({ a <- rt_analysis(d, list(covs = "X001")); a$run() }, "implementation limitation")
  testthat::expect_error({ a <- rt_analysis(d, list(cvNfolds = 100)); a$run() }, "folds")
  testthat::local_mocked_bindings(ridgetrace = function(...) stop("test engine failure"), .package = "GCEstim")
  testthat::expect_error({ a <- rt_analysis(d); a$run() }, "Ridge estimation failed.*test engine failure")
})

testthat::test_that("ridge core plots render and persisted cache never refits for presentation", {
  a <- rt_analysis(rt_data(), list(plotSample = "show")); a$run()
  state <- a$results$fitCache$state
  testthat::expect_true(rt_plain(state))
  testthat::local_mocked_bindings(ridgetrace = function(...) stop("Unexpected ridge refit"), .package = "GCEstim")
  for (changes in list(list(selection = "1se"), list(selection = "manual", manualIndex = 3L),
      list(axisScale = "linear", extrema = TRUE, facet = TRUE), list(showIntercept = FALSE),
      list(plotCV = FALSE), list(tableRows = 2L))) {
    replay <- rt_analysis(rt_data(), changes)
    replay$results$fitCache$.__enclos_env__$super$fromProtoBuf(state_message(state), character(), character())
    replay$run()
    testthat::expect_identical(replay$results$fitCache$state, state)
    testthat::expect_match(replay$results[[".fitTime"]]$content, "Previous fit reused")
    testthat::expect_equal(nrow(replay$results$lambdaResults$asDF), min(if (is.null(changes$tableRows)) 100 else changes$tableRows, 8))
  }
  device <- tempfile(fileext = ".pdf"); grDevices::pdf(device)
  tryCatch(for (key in c("plotPaths", "plotCV", "plotSample")) {
    testthat::expect_true(a$.__enclos_env__$private[[paste0(".",key)]](a$results[[key]],
      ggplot2::theme_minimal(), list(color = c("blue","green"))))
  }, finally = { grDevices::dev.off(); unlink(device) })
  before <- a$results$coefficients$asDF; a$run()
  testthat::expect_equal(a$results$coefficients$asDF, before)
  clear <- a$results$fitCache$.__enclos_env__$private$.clearWith
  testthat::expect_true(all(c("data", "standardize", "cv", "lambdaValues", "seed", "errorMeasure") %in% clear))
  testthat::expect_false(any(c("selection", "manualIndex", "plotPaths", "axisScale", "displayTerms", "tableRows") %in% clear))
})


testthat::test_that("all ridge plots render from persisted numerical state", {
  keys <- c("plotPaths", "plotCV", "plotSample", "plotRanges", "plotSensitivity", "plotFolds", "plotStability")
  d <- rt_data(); d$X005 <- factor(d$X005)
  a <- rt_analysis(d, c(list(factors = "X005"), stats::setNames(lapply(keys, function(k) if (k == "plotSample") "show" else TRUE), keys)))
  a$run(); cache <- a$results$fitCache$state
  testthat::expect_false(".fitWarning" %in% names(a$results))
  for (key in keys) {
    testthat::expect_true(a$results[[key]]$visible, info = key)
    testthat::expect_true(rt_plain(a$results[[key]]$state))
    testthat::expect_null(a$results[[key]]$state$data)
    testthat::expect_null(a$results[[key]]$state$selection$mean)
  }
  testthat::local_mocked_bindings(ridgetrace = function(...) stop("Unexpected refit"),
    plot.ridgetrace = function(...) stop("Plot extraction forbidden"), .package = "GCEstim")
  device <- tempfile(fileext = ".pdf"); grDevices::pdf(device, width = 9, height = 7)
  tryCatch(for (axis in c("log", "linear")) {
    replay <- rt_analysis(d, c(list(factors = "X005", axisScale = axis,
      extrema = TRUE, facet = TRUE, showIntercept = TRUE), stats::setNames(lapply(keys, function(k) if (k == "plotSample") "show" else TRUE), keys)))
    replay$results$fitCache$.__enclos_env__$super$fromProtoBuf(state_message(cache), character(), character())
    replay$run()
    testthat::expect_false(".fitWarning" %in% names(replay$results))
    for (key in keys) {
      testthat::expect_true(replay$.__enclos_env__$private[[paste0(".", key)]](
        replay$results[[key]], ggplot2::theme_minimal(), list(color = c("blue","green","red"))), info = key)
    }
  }, finally = { grDevices::dev.off(); unlink(device) })
  s <- a$results$plotSensitivity$state
  p <- GCEstimjmv:::.gceRTPlot(cache, s, NULL, list())
  b <- cache$fit$coefficients
  expected <- (b[, -1, drop = FALSE] - b[, -ncol(b), drop = FALSE]) / rep(diff(log10(cache$fit$lambda)), each = nrow(b))
  testthat::expect_equal(p$data$estimate, as.numeric(expected))
  invalid <- rt_analysis(d, list(factors = "X005", displayTerms = "unknown", plotSample = "show"))
  invalid$results$fitCache$setState(cache); invalid$run()
  testthat::expect_false(invalid$results$plotPaths$visible)
  testthat::expect_true(invalid$results$plotCV$visible)
  testthat::expect_true(invalid$results$plotSample$visible)
  testthat::expect_match(invalid$results[[".fitWarning"]]$content, "Unknown coefficient")
  testthat::expect_equal(invalid$results$coefficients$asDF, a$results$coefficients$asDF)
})

testthat::test_that("ridge nonfinite returned errors preserve usable coefficient results", {
  direct <- GCEstim::ridgetrace(y ~ X001 + X002, rt_data(), lambda.n = 8, cv.nfolds = 3)
  direct$error.lambda[] <- Inf; direct$error.lambda.cv[,] <- NaN
  testthat::local_mocked_bindings(ridgetrace = function(...) direct, .package = "GCEstim")
  a <- rt_analysis(rt_data()); a$run()
  testthat::expect_equal(nrow(a$results$coefficients$asDF), 3L)
  testthat::expect_true(all(is.na(a$results$coefficients$asDF$estimate)))
  testthat::expect_true(a$results$plotPaths$visible)
  testthat::expect_match(a$results[[".fitWarning"]]$content, "nonfinite")
  manual <- rt_analysis(rt_data(), list(selection = "manual", manualIndex = 2L)); manual$run()
  testthat::expect_equal(manual$results$coefficients$asDF$estimate, as.numeric(direct$coef.lambda[,2]))
})

testthat::test_that("ridge UI and dependencies retain the approved scope", {
  path <- function(x) normalizePath(file.path(testthat::test_path(), "..", "..", x), winslash = "/")
  a <- readLines(path("jamovi/ridgetrace.a.yaml"))
  testthat::expect_false(any(grepl("^      - (MAPE|MASE)$", a)))
  analysis <- rt_analysis(rt_data()); deps <- analysis$results$fitCache$.__enclos_env__$private$.clearWith
  for (name in c("data", "run", "dep", "covs", "factors", "intercept", "standardize", "penalizeIntercept", "lambdaN", "lambdaMin", "lambdaMax", "customGrid", "lambdaValues", "cv", "cvNfolds", "seed", "errorMeasure")) {
    testthat::expect_true(name %in% deps, info = name)
    results <- rt_analysis(rt_data())$results
    results$fitCache$.__enclos_env__$super$fromProtoBuf(state_message(list(version = 1L)), name, character())
    testthat::expect_null(results$fitCache$state, info = name)
  }
  for (name in c("selection", "manualIndex", "showLambdaTable", "tableRows", "axisScale", "displayTerms", "showIntercept", "extrema", "facet", "seBars", "foldMean", "plotPaths", "plotStability"))
    testthat::expect_false(name %in% deps, info = name)
})


testthat::test_that("ridge supported error measures are the engine's unchanged errors", {
  d <- rt_data()
  for (measure in c("RMSE", "MSE", "MAE", "sMAPE")) {
    a <- rt_analysis(d, list(errorMeasure = measure)); a$run()
    direct <- GCEstim::ridgetrace(y ~ X001 + X002, d, lambda.n = 8, cv.nfolds = 3, errormeasure = measure)
    testthat::expect_equal(a$results$fitCache$state$fit$sample, direct$error.lambda)
    testthat::expect_equal(unname(a$results$fitCache$state$fit$cvErrors), unname(as.matrix(direct$error.lambda.cv)))
  }
})

testthat::test_that("ridge names, single-predictor no-CV and in-sample visibility are preserved", {
  d <- rt_data(); names(d)[names(d) == "X001"] <- "Signal, X001"
  a <- rt_analysis(d, list(covs = c("Signal, X001", "X002"), cv = FALSE)); a$run()
  testthat::expect_true(any(grepl("Signal, X001", a$results$coefficients$asDF$term, fixed = TRUE)))
  f <- a$results$fitCache$state$fit
  term <- rownames(f$coefficients)[2]
  quoted <- paste0('"', term, '"')
  testthat::expect_identical(GCEstimjmv:::.gceRTTerms(f, FALSE, quoted), term)
  single <- rt_analysis(rt_data(), list(covs = "X001", cv = FALSE)); single$run()
  testthat::expect_equal(nrow(single$results$coefficients$asDF), 2)
  hidden <- rt_analysis(rt_data(), list(cv = FALSE, plotSample = "hide")); hidden$run()
  testthat::expect_false(hidden$results$plotSample$visible)
})


testthat::test_that("ridge reference styles and extrema match TARW conventions", {
  a <- rt_analysis(rt_data(), list(extrema = TRUE, plotRanges = TRUE)); a$run()
  cache <- a$results$fitCache$state
  p <- GCEstimjmv:::.gceRTPlot(cache, a$results$plotCV$state, NULL, list())
  vertical <- which(vapply(p$layers, function(x) inherits(x$geom, "GeomVline"), logical(1)))
  testthat::expect_length(vertical, 1L)
  layer <- ggplot2::ggplot_build(p)$data[[vertical]]
  testthat::expect_identical(layer$linetype, c("dotted", "dashed"))
  testthat::expect_true(all(layer$colour == "black"))
  testthat::expect_equal(layer$linewidth, rep(.5, 2))
  s <- a$results$plotCV$state$selection
  testthat::expect_equal(layer$xintercept, log10(cache$fit$lambda[c(s$minimum, s$oneSE)]))
  testthat::expect_identical(p$scales$get_scales("linetype")$breaks, c("Minimum error", "1-SE rule"))
  p <- GCEstimjmv:::.gceRTPlot(cache, a$results$plotPaths$state, NULL, list(color = c("blue", "green")))
  points <- Filter(function(x) inherits(x$geom, "GeomPoint"), p$layers)
  testthat::expect_length(points, 2L)
  testthat::expect_equal(unname(vapply(points, function(x) x$aes_params$shape, numeric(1))), c(24,25))
  testthat::expect_equal(unname(vapply(points, function(x) x$aes_params$size, numeric(1))), c(1,1))
  for (j in seq_along(points)) {
    b <- cache$fit$coefficients
    indices <- apply(b, 1L, if (j == 1L) which.max else which.min)
    testthat::expect_equal(points[[j]]$data$lambda, cache$fit$lambda[indices])
    testthat::expect_equal(points[[j]]$data$estimate, as.numeric(b[cbind(seq_len(nrow(b)), indices)]))
    testthat::expect_identical(as.character(points[[j]]$data$term), rownames(b))
  }
  built <- ggplot2::ggplot_build(p)
  pointIndices <- which(vapply(p$layers, function(x) inherits(x$geom, "GeomPoint"), logical(1)))
  colours <- GCEstimjmv:::.gceTermColours(rownames(cache$fit$coefficients), list(color = c("blue","green")))
  for (i in pointIndices) testthat::expect_equal(built$data[[i]]$colour, unname(colours[rownames(cache$fit$coefficients)]))
  ranges <- GCEstimjmv:::.gceRTPlot(cache, a$results$plotRanges$state, NULL, list())
  testthat::expect_equal(ranges$data$estimate, as.numeric(cache$fit$coefficients[, s$index]))
  testthat::expect_match(ranges$labels$caption, "Points: selected[[:space:]]+coefficients")
})

testthat::test_that("ridge registers seven plots with display defaults and no retired controls", {
  a <- rt_analysis(rt_data()); a$run()
  keys <- c("plotPaths", "plotRanges", "plotSensitivity", "plotCV", "plotSample", "plotFolds", "plotStability")
  retired <- c("plotHeatmap", "plotSelected", "normalize")
  testthat::expect_true(a$options$showIntercept)
  testthat::expect_equal(a$options$seed, 230676L)
  testthat::expect_true(a$options$intercept)
  testthat::expect_false(any(retired %in% names(a$options$options)))
  testthat::expect_false(any(retired %in% names(a$results)))
  for (key in keys) testthat::expect_true(inherits(a$results[[key]], "Image"), info = key)
  base <- file.path(testthat::test_path(), "..", "..")
  definitions <- yaml::read_yaml(file.path(base, "jamovi/ridgetrace.r.yaml"))
  images <- Filter(function(x) x$type == "Image", definitions$items)
  testthat::expect_setequal(vapply(images, function(x) x$name, character(1)), keys)
  testthat::expect_length(images, 7L)
  ui <- yaml::read_yaml(file.path(base, "jamovi/ridgetrace.u.yaml"))
  flatten <- function(x) c(list(x), unlist(lapply(x$children, flatten), recursive = FALSE))
  controls <- flatten(ui)
  seed <- Filter(function(x) identical(x$name, "seed"), controls)[[1L]]
  testthat::expect_identical(seed$width, "largest")
  testthat::expect_false(any(vapply(controls, function(x) (!is.null(x$name) && x$name %in% retired) || (!is.null(x$enable) && grepl("plotHeatmap", x$enable)), logical(1))))
  testthat::expect_false(grepl("derived by GCEstimjmv|1-SE reference is derived", a$results[[".method"]]$content))
  for (key in c("plotPaths", "plotRanges", "plotSensitivity", "plotStability")) {
    o <- a$results$plotPaths$state
    o$key <- key
    if (key == "plotStability") o$options$showIntercept <- TRUE
    p <- GCEstimjmv:::.gceRTPlot(a$results$fitCache$state, o, NULL, list())
    testthat::expect_true("(Intercept)" %in% as.character(p$data$term), info = key)
  }
})

testthat::test_that("saved obsolete ridge options and images are retired without refitting", {
  a <- rt_analysis(rt_data(), list(plotRanges = TRUE)); a$run()
  cache <- a$results$fitCache$state
  legacyImage <- a$results$plotPaths$state
  legacyImage$options$normalize <- TRUE
  pb <- list(state = raw(), group = list(elements = list(
    c(list(name = "fitCache"), state_message(cache)),
    list(name = "plotHeatmap", state = raw()), list(name = "plotSelected", state = raw()))))
  replay <- rt_analysis(rt_data(), list(plotPaths = FALSE, selection = "1se"))
  replay$results$fromProtoBuf(pb, character(), character())
  replay$results$plotPaths$.__enclos_env__$super$fromProtoBuf(state_message(legacyImage), character(), character())
  testthat::expect_identical(replay$results$fitCache$state, cache)
  retired <- c("plotHeatmap", "plotSelected", "normalize")
  testthat::local_mocked_bindings(parseOptionPB = function(pb) pb$value, .package = "jmvcore")
  enumName <- "jamovi.coms.AnalysisOption.Other"
  hadEnum <- exists(enumName, envir = globalenv(), inherits = FALSE)
  if (hadEnum) oldEnum <- get(enumName, envir = globalenv())
  withr::defer(if (hadEnum) assign(enumName, oldEnum, envir = globalenv()) else rm(list = enumName, envir = globalenv()))
  assign(enumName, list(NONE = 2L), envir = globalenv())
  replay$options$fromProtoBuf(list(names = retired, options = lapply(retired, function(x) list(value = TRUE, o = 1L))))
  testthat::expect_true(all(vapply(retired, replay$options$has, logical(1))))
  testthat::local_mocked_bindings(ridgetrace = function(...) stop("Unexpected refit"), .package = "GCEstim")
  replay$run()
  testthat::expect_false(any(vapply(retired, replay$options$has, logical(1))))
  testthat::expect_true(all(vapply(replay$options$asProtoBuf()$options, function(x) x$o == 2L, logical(1))))
  testthat::expect_false(any(retired %in% names(replay$results)))
  testthat::expect_null(replay$results$plotPaths$state$options$normalize)
  testthat::expect_identical(replay$results$fitCache$state, cache)
  testthat::expect_match(replay$results[[".fitTime"]]$content, "Previous fit reused")
  testthat::expect_equal(replay$results$coefficientRanges$asDF$selected,
    as.numeric(cache$fit$coefficients[, GCEstimjmv:::.gceRTSelection(cache$fit, "1se", 1, TRUE)$index]))
  for (show in c(FALSE, TRUE)) {
    fresh <- rt_analysis(rt_data(), list(showIntercept = show, selection = "manual", manualIndex = 2L))
    fresh$results$fitCache$setState(cache); fresh$run()
    testthat::expect_identical(fresh$results$fitCache$state, cache)
    testthat::expect_equal(fresh$results$coefficients$asDF$estimate, as.numeric(cache$fit$coefficients[,2]))
    testthat::expect_false(any(retired %in% names(fresh$results$plotPaths$state$options)))
  }
})
