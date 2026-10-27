# Optional overrides: GCESTIMJMV_R_LIBRARY and JAMOVI_HOME.
r_series <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][1L], sep = ".")
r_library <- Sys.getenv("GCESTIMJMV_R_LIBRARY", file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", r_series))
.libPaths(c(r_library, .libPaths()))
jamovi_candidates <- list.dirs(Sys.getenv("ProgramFiles", "C:/Program Files"), recursive = FALSE)
jamovi_candidates <- jamovi_candidates[grepl("^jamovi [0-9]", basename(jamovi_candidates))]
options(jamovi_home = Sys.getenv("JAMOVI_HOME", if (length(jamovi_candidates)) tail(sort(jamovi_candidates), 1L) else ""))
args <- commandArgs(TRUE)
if ("prepare" %in% args) {
  status <- jmvtools::prepare()
  if (!identical(status, 0L)) stop("jamovi prepare failed: ", status)
}
if ("core" %in% args) {
  devtools::load_all()
  d <- read.csv("tests/testthat/_fixtures/dataThesis.csv")
  a <- RidgeTrace(data=d, dep=y, covs=vars(X001,X002), factors=NULL, run=TRUE, lambdaN=8, cvNfolds=3, showLambdaTable=TRUE)
  print(a$modelSummary$asDF); print(a$coefficients$asDF); print(a$lambdaResults$asDF)
}
if ("test" %in% args) {
  dir.create("build/ridge-validation", showWarnings = FALSE, recursive = TRUE)
  sink("build/ridge-validation/tests.log", split = TRUE)
  result <- devtools::test(stop_on_failure = TRUE)
  sink()
  frame <- as.data.frame(result)
  frame$result <- NULL
  utils::write.csv(frame, "build/ridge-validation/test-results.csv", row.names = FALSE)
}
if ("ridge" %in% args) devtools::test(filter = "ridgetrace", stop_on_failure = TRUE)

if ("plots" %in% args) {
  devtools::load_all()
  d <- read.csv("tests/testthat/_fixtures/dataThesis.csv")
  options <- RidgeTraceOptions$new(dep = "y", covs = c("X001", "X002"), factors = "X005", run = TRUE,
    lambdaN = 20, cvNfolds = 3, plotSample = "show", plotRanges = TRUE,
    plotSensitivity = TRUE, plotFolds = TRUE, plotStability = TRUE)
  analysis <- RidgeTraceClass$new(options = options, data = d); analysis$run()
  dir.create("build/ridge-validation", showWarnings = FALSE)
  for (key in c("plotPaths", "plotCV", "plotSample", "plotRanges", "plotSensitivity", "plotFolds", "plotStability")) {
    image <- analysis$results[[key]]
    p <- GCEstimjmv:::.gceRTPlot(analysis$results$fitCache$state, image$state,
      ggplot2::theme_minimal(), list(color = c("#346da4", "#e1812c", "#399c61", "#b95696")))
    ggplot2::ggsave(file.path("build/ridge-validation", paste0(key, ".png")), p, width = 8, height = 6, dpi = 120)
  }
}

if ("build" %in% args) {
  build_ridge <- function() {
    dir.create("build/ridge-validation", showWarnings = FALSE, recursive = TRUE)
    output <- normalizePath("build/ridge-validation", winslash = "/", mustWork = TRUE)
    backup <- file.path(output, "temp-backup")
    dir.create(backup, showWarnings = FALSE)
    originals <- list.files("temp", all.files = TRUE, no.. = TRUE, full.names = TRUE)
    if (length(originals) && !all(file.copy(originals, backup, recursive = TRUE, copy.date = TRUE))) stop("Could not preserve compiler scratch files")
    before <- tools::md5sum(originals)
    on.exit({
      saved <- file.path(backup, basename(originals))
      dir.create("temp", showWarnings = FALSE)
      if (length(saved) && !all(file.copy(saved, "temp", recursive = TRUE, overwrite = TRUE, copy.date = TRUE))) stop("Scratch restoration failed")
      if (!identical(unname(tools::md5sum(originals)), unname(before))) stop("Scratch file hashes changed")
    })
    jmc <- file.path(system.file("node_modules/jamovi-compiler", package = "jmvtools"), "index.js")
    status <- system2(Sys.which("node"), c(shQuote(jmc), "--build", shQuote(normalizePath(".", winslash="/")),
      "--home", shQuote(getOption("jamovi_home")), "--skip-remotes", "--skip-deps",
      "--jmo", shQuote(file.path(output, "GCEstimjmv-ridge-development.jmo"))))
    if (!identical(status, 0L)) stop("Development build failed: ", status)
  }
  build_ridge()
}

if ("bundled" %in% args) {
  if (isNamespaceLoaded("GCEstimjmv")) stop("Run bundled in a fresh Rscript invocation.")
  .libPaths(c("build/R4.5.0-x64-win64", file.path(getOption("jamovi_home"), "Resources/modules/base/R"),
    file.path(getOption("jamovi_home"), "Resources/modules/jmv/R")))
  library(GCEstimjmv)
  stopifnot(normalizePath(getNamespaceInfo(asNamespace("GCEstimjmv"), "path")) == normalizePath("build/R4.5.0-x64-win64/GCEstimjmv"))
  d <- read.csv("tests/testthat/_fixtures/dataThesis.csv")
  a <- GCEstimjmv::RidgeTrace(data=d, dep=y, covs=vars(X001,X002), factors=NULL, run=TRUE, lambdaN=8, cvNfolds=3)
  stopifnot(nrow(a$coefficients$asDF) == 3L, a$plotPaths$visible, a$plotCV$visible)
  cat("Bundled installed package smoke check passed\n")
}

if ("archive" %in% args) {
  archive <- "build/ridge-validation/GCEstimjmv-ridge-development.jmo"
  members <- utils::unzip(archive, list = TRUE)$Name
  read_member <- function(name) {
    connection <- unz(archive, name)
    on.exit(close(connection))
    readLines(connection, warn = FALSE)
  }
  meta <- read_member("GCEstimjmv/jamovi.yaml")
  description <- read_member("GCEstimjmv/R/GCEstimjmv/DESCRIPTION")
  stopifnot(any(grepl("name: RidgeTrace", meta, fixed = TRUE)),
    any(grepl("/ridgetrace.js$", members)), any(grepl("utils", description, fixed = TRUE)))
  cat("Development JMO menu, compiled UI and package metadata verified\n")
}
