# Local preparation only. No commit, publication or estimator changes.
# Reuse the validated workstation/library discovery and test/build conventions.
source("tools/validate-ridge.R")
output <- "build/local-1.1.0"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
expected <- c("WGCE", "TARWGCE", "TASWGCE", "RidgeTrace")
archive <- normalizePath(".", winslash = "/")
archive <- file.path(archive, "GCEstimjmv_1.1.0.jmo")
read_member <- function(name) {
  connection <- unz(archive, paste0("GCEstimjmv/", name))
  on.exit(close(connection))
  readLines(connection, warn = FALSE)
}
verify_archive <- function() {
  meta <- yaml::yaml.load(paste(read_member("jamovi.yaml"), collapse = "\n"))
  full <- yaml::yaml.load(paste(read_member("jamovi-full.yaml"), collapse = "\n"))
  connection <- textConnection(read_member("R/GCEstimjmv/DESCRIPTION"))
  on.exit(close(connection))
  description <- read.dcf(connection)
  stopifnot(as.character(meta$version) == "1.1.0", as.character(full$version) == "1.1.0",
    description[1L, "Version"] == "1.1.0",
    setequal(vapply(meta$analyses, function(x) x$name, character(1)), expected),
    setequal(vapply(full$analyses, function(x) x$name, character(1)), expected))
  members <- utils::unzip(archive, list = TRUE)$Name
  stopifnot(all(paste0("GCEstimjmv/ui/", c("wgce", "tarwgce", "taswgce", "ridgetrace"), ".js") %in% members))
  ridge <- full$analyses[[which(vapply(full$analyses, function(x) x$name == "RidgeTrace", logical(1)))]]
  options <- stats::setNames(ridge$options, vapply(ridge$options, function(x) x$name, character(1)))
  stopifnot(!any(c("plotHeatmap", "plotSelected", "normalize") %in% names(options)),
    isTRUE(options$showIntercept$default), setequal(unlist(options$errorMeasure$options), c("RMSE", "MSE", "MAE", "sMAPE")))
  writeLines(c("Embedded jamovi.yaml: 1.1.0", "Embedded jamovi-full.yaml: 1.1.0",
    "Embedded R DESCRIPTION: 1.1.0", paste("Analyses:", paste(expected, collapse = ", ")),
    "All four compiled UI definitions present", "Four approved error measures verified"), file.path(output, "archive-verification.txt"))
  cat("Embedded package metadata, analyses and compiled UI verified\n")
}
if ("inspect110" %in% args)
  cat(readLines(file.path(system.file("node_modules/jamovi-compiler", package = "jmvtools"), "installer.js")), sep = "\n")
if ("preflight110" %in% args) {
  description <- read.dcf("DESCRIPTION")
  meta <- yaml::read_yaml("jamovi/0000.yaml")
  stopifnot(description[1L,"Version"] == "1.1.0", as.character(meta$version) == "1.1.0",
    setequal(vapply(meta$analyses, function(x) x$name, character(1)), expected))
  definitions <- list.files("jamovi", pattern = "[.]yaml$", full.names = TRUE)
  invisible(lapply(definitions, yaml::read_yaml))
  imports <- trimws(gsub(" *[(].*[)]", "", strsplit(description[1L,"Imports"], ",")[[1L]]))
  stopifnot(setequal(imports, c("jmvcore", "R6", "GCEstim", "stats", "ggplot2", "grDevices", "utils")),
    grepl("GCEstim (>= 1.1.0)", description[1L,"Imports"], fixed = TRUE),
    description[1L,"Remotes"] == "jorgevazcabral/GCEstim@v1.1.0")
  a <- yaml::read_yaml("jamovi/ridgetrace.a.yaml")
  r <- yaml::read_yaml("jamovi/ridgetrace.r.yaml")
  options <- stats::setNames(a$options, vapply(a$options, function(x) x$name, character(1)))
  images <- Filter(function(x) x$type == "Image", r$items)
  stopifnot(length(images) == 7L,
    setequal(unlist(options$errorMeasure$options), c("RMSE", "MSE", "MAE", "sMAPE")))
  backend <- paste(readLines("R/ridgetrace.b.R"), collapse = "\n")
  stopifnot(grepl("implementation limitation, not a statistical restriction", backend, fixed = TRUE),
    grepl("complete dataset before constructing CV folds", backend, fixed = TRUE),
    grepl("GCEstim::ridgetrace(", backend, fixed = TRUE))
  cache <- Filter(function(x) x$name == "fitCache", r$items)[[1L]]
  stopifnot(!any(c("selection", "manualIndex", "axisScale", "plotPaths", "tableRows") %in% unlist(cache$clearWith)))
  writeLines(c(paste("Analyses:", paste(expected, collapse=", ")),
    paste("YAML files parsed:", length(definitions)), "Ridge Trace images: 7",
    "Dependencies unchanged; GCEstim pinned to v1.1.0", "Guards, CV note and cache dependencies verified"), file.path(output,"preflight.txt"))
  cat("Version 1.1.0 preflight passed\n")
}
if ("tests110" %in% args) {
  status <- system2(file.path(R.home("bin"), "Rscript.exe"), c("--vanilla", "tools/validate-ridge.R", "test"))
  if (!identical(status, 0L)) stop("Complete suite failed: ", status)
  stopifnot(file.copy("build/ridge-validation/tests.log", output, overwrite = TRUE),
    file.copy("build/ridge-validation/test-results.csv", output, overwrite = TRUE))
  results <- utils::read.csv(file.path(output,"test-results.csv"))
  stopifnot(sum(results$passed) >= 1879L, sum(results$failed) == 0L,
    sum(results$warning) == 0L, !any(results$error), !any(results$skipped))
  cat("Complete suite:", sum(results$passed), "passed; Ridge Trace:", sum(results$passed[results$file == "test-ridgetrace.R"]), "; zero failures, test warnings or skips\n")
}
if ("build110" %in% args) {
  build_local <- function() {
    backup <- file.path(output, "temp-backup")
    dir.create(backup, showWarnings = FALSE)
    originals <- list.files("temp", all.files = TRUE, no.. = TRUE, full.names = TRUE)
    if (length(originals) && !all(file.copy(originals, backup, recursive = TRUE, overwrite = TRUE, copy.date = TRUE))) stop("Could not preserve scratch files")
    before <- tools::md5sum(originals)
    on.exit({
      dir.create("temp", showWarnings = FALSE)
      saved <- file.path(backup, basename(originals))
      if (length(saved) && !all(file.copy(saved, "temp", recursive = TRUE, overwrite = TRUE, copy.date = TRUE))) stop("Scratch restoration failed")
      stopifnot(identical(unname(tools::md5sum(originals)), unname(before)))
    })
    jmc <- file.path(system.file("node_modules/jamovi-compiler", package = "jmvtools"), "index.js")
    status <- system2(Sys.which("node"), c(shQuote(jmc), "--build", shQuote(normalizePath(".", winslash="/")),
      "--home", shQuote(getOption("jamovi_home")), "--skip-remotes", "--skip-deps", "--jmo", shQuote(archive)))
    if (!identical(status, 0L)) stop("Local build failed: ", status)
  }
  build_local()
  verify_archive()
}
if ("archive110" %in% args) verify_archive()
if ("smoke110" %in% args) {
  if (isNamespaceLoaded("GCEstimjmv")) stop("Smoke test requires a fresh R process")
  .libPaths(c("build/R4.5.0-x64-win64", file.path(getOption("jamovi_home"), "Resources/modules/base/R"),
    file.path(getOption("jamovi_home"), "Resources/modules/jmv/R")))
  library(GCEstimjmv)
  stopifnot(normalizePath(getNamespaceInfo(asNamespace("GCEstimjmv"), "path")) == normalizePath("build/R4.5.0-x64-win64/GCEstimjmv"),
    as.character(utils::packageVersion("GCEstimjmv")) == "1.1.0",
    as.character(utils::packageVersion("GCEstim")) == "1.1.0")
  stopifnot(all(expected %in% getNamespaceExports("GCEstimjmv")))
  d <- utils::read.csv("tests/testthat/_fixtures/dataThesis.csv")
  a <- GCEstimjmv::RidgeTrace(data = d, dep = y, covs = vars(X001, X002), factors = NULL,
    run = TRUE, lambdaN = 8, cvNfolds = 3)
  direct <- GCEstim::ridgetrace(y ~ X001 + X002, d, lambda.n = 8, cv.nfolds = 3)
  selected <- which.min(colMeans(direct$error.lambda.cv))
  stopifnot(isTRUE(all.equal(a$coefficients$asDF$estimate, as.numeric(direct$coef.lambda[,selected]))),
    a$plotPaths$visible, a$plotCV$visible)
  writeLines(c("Fresh-process compiled GCEstimjmv 1.1.0 loaded", "GCEstim 1.1.0 available",
    "All four analysis entry points exported", "Ridge Trace example matched direct engine coefficients"), file.path(output, "smoke-check.txt"))
  cat("Compiled v1.1.0 fresh-process loading and numerical smoke check passed\n")
}
if ("install110" %in% args) {
  verify_archive()
  installer <- file.path(system.file("node_modules/jamovi-compiler", package = "jmvtools"), "installer.js")
  installer <- paste0("file:///", normalizePath(installer, winslash = "/"))
  helper <- file.path(output, "install-local.mjs")
  writeLines(c(paste0("import installer from ", jsonlite::toJSON(installer, auto_unbox=TRUE), ";"),
    paste0("installer.install(", jsonlite::toJSON(archive, auto_unbox=TRUE), ", ",
      jsonlite::toJSON(getOption("jamovi_home"), auto_unbox=TRUE), ");")), helper)
  status <- system2(Sys.which("node"), shQuote(helper), stdout = file.path(output,"installation.log"), stderr = file.path(output,"installation.log"))
  cat(readLines(file.path(output,"installation.log")), sep="\n")
  if (!identical(status, 0L)) stop("Supported local installer failed: ", status)
  cat("\nSupported jamovi local installer returned success\n")
}
