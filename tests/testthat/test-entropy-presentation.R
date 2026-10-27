
entropy_plot <- function(state) {
  fun <- GCEstimjmv:::.gceRenderDiagnostic
  env <- new.env(parent = environment(fun))
  env$print <- function(x) { env$p <- x; invisible(x) }
  environment(fun) <- env
  g <- jmvcore:::getGlobalTheme("default", "jmv")
  expect_true(suppressMessages(fun(list(state = state), g$ggtheme, g$theme)))
  list(plot = env$p, built = ggplot2::ggplot_build(env$p))
}

entropy_data <- function() {
  d <- read.csv(test_path("_fixtures", "dataThesis.csv"))
  names(d)[c(1,2,5)] <- c("Z signal", "A signal", "F group")
  d[["F group"]] <- factor(d[["F group"]]); d
}

entropy_analysis <- function(kind) {
  args <- list(dep = "y", covs = c("Z signal", "A signal"), factors = "F group", run = TRUE,
    M = "3", J = "3", weight = "0.5", cvNfolds = 5L, supportSignalVectorN = 5L,
    twostepsN = 0L, plot1 = TRUE, plot3 = TRUE, plot4 = TRUE)
  if (kind == "TARWGCE") args$ridgeLambdaN <- 20L
  ns <- asNamespace("GCEstimjmv")
  get(paste0(kind, "Class"), ns)$new(
    options = do.call(get(paste0(kind, "Options"), ns)$new, args), data = entropy_data())
}

test_that("failed first and intermediate supports cannot discard finite entropy geometry", {
  for (failed in c(1,3)) for (kind in c("TARWGCE", "TASWGCE")) {
    ids <- c("Z signal", "A signal", "\u0060F group\u0060B")
    d <- expand.grid(predictor = ids, support = c(.3,.8572,2.4495,6.9993,20),
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    d$nepk <- seq_len(nrow(d)) / (nrow(d) + 1)
    d$estimate <- NA_real_; d$unused <- NaN
    d$nepk[d$support == unique(d$support)[failed]] <- NaN
    good <- d[is.finite(d$nepk), ]
    refs <- list(min = 2.4495, oneSE = .8572, manual = .6)
    source <- list(which = 3L, data = d[sample(seq_len(nrow(d))), ], supportLines = refs)
    state <- list(which = 4L, coefficientTerms = ids,
      data = good[good$support != good$support[1], ], supportLines = refs, overlays = list())
    before <- serialize(list(source, state), NULL)
    recovered <- GCEstimjmv:::.gceEntropyState(state, source)
    expect_true(recovered$entropyComplete)
    expect_null(recovered$entropyNotice)
    expect_equal(nrow(recovered$data), 12L)
    expect_false(anyDuplicated(recovered$data[c("support","predictor")]) > 0L)
    aligned <- merge(good, recovered$data, by = c("support","predictor"))
    expect_identical(aligned$nepk.x, aligned$nepk.y)
    expect_true(all(is.na(recovered$data$estimate)))
    expect_identical(serialize(list(source,state),NULL), before)
    p <- entropy_plot(recovered)
    expect_equal(nrow(p$built$data[[1]]), 12L)
    refsBuilt <- Filter(function(x) "xintercept" %in% names(x), p$built$data)
    expect_identical(vapply(refsBuilt,function(x)x$xintercept[1],numeric(1)),c(2.4495,.8572,.6))
    expect_identical(vapply(refsBuilt,function(x)x$linetype[1],character(1)),c("dotted","dashed","dotted"))
    expect_identical(vapply(refsBuilt,function(x)x$colour[1],character(1)),c("black","black","green4"))
  }
})

test_that("legacy entropy recovery rejects mismatched identities, stages and coordinates", {
  d <- data.frame(support = rep(c(1,2),each=2), predictor = rep(c("Z","A"),2),
                  estimate = NA_real_, nepk = c(.2,.3,.6,.7))
  s <- list(which=4L, coefficientTerms=c("Z","A"),data=d[3:4,],overlays=list(),
            supportLines=list(min=2))
  source <- list(which=3L,data=d,supportLines=s$supportLines)
  for (bad in c("identifier","coordinate","entropy","duplicate","stage","reference")) {
    x <- source
    if (bad=="identifier") x$data$predictor[1] <- "other"
    if (bad=="coordinate") x$data$support[3] <- 3
    if (bad=="entropy") x$data$nepk[3] <- .8
    if (bad=="duplicate") x$data <- rbind(x$data,x$data[1,])
    if (bad=="stage") x$which <- 6L
    if (bad=="reference") x$supportLines$min <- 1
    result <- GCEstimjmv:::.gceEntropyState(s,x)
    expect_false(isTRUE(result$entropyComplete))
    expect_identical(result$data,s$data)
    expect_match(result$entropyNotice,"Only available observations")
    p <- entropy_plot(result)
    expect_match(p$plot$labels$caption,"missing values")
  }
  legacy <- unserialize(serialize(s,NULL))
  image <- GCEstimjmv:::.gcePlotImage(list(state=legacy),list(plots=list(plot3=source)))
  expect_equal(nrow(image$state$data),4L)
  missing <- GCEstimjmv:::.gcePlotImage(list(state=legacy),list())
  expect_equal(nrow(missing$state$data),2L)
  expect_match(missing$state$entropyNotice,"unavailable")
})

for (kind in c("TARWGCE","TASWGCE")) local({
  kind <- kind
  test_that(paste(kind,"preserves fitted entropy and repairs serialized caches without estimation"),{
    full <- NULL; original <- GCEstim::cv.lmgce
    local_mocked_bindings(cv.lmgce=function(...) {full <<- original(...);full},.package="GCEstim")
    a <- entropy_analysis(kind);a$run()
    cache <- a$results$fitCache$state
    state <- GCEstimjmv:::.gcePlotImage(a$results$plot4,cache)$state
    n <- if(kind=="TARWGCE")20L else 25L
    expect_equal(nrow(state$data),n)
    expect_true(state$entropyComplete)
    best <- full$best;ids <- names(coef(best))
    expected <- do.call(rbind,lapply(colnames(best$results$nocv$support.results$coef.matrix.ss),function(s) {
      q <- best$results$nocv$support.results[[s]]
      data.frame(predictor=ids,support=as.numeric(s),nepk=unname(q$nepk[ids]))
    }))
    aligned <- merge(state$data,expected,by=c("support","predictor"))
    expect_equal(aligned$nepk.x,aligned$nepk.y,tolerance=1e-12)
    expect_equal(nrow(aligned),n)
    point <- state$overlays[[1]]
    vertices <- merge(point,state$data,by.x=c("predictor","x"),by.y=c("predictor","support"))
    expect_equal(nrow(vertices),length(ids))
    expect_equal(vertices$y,vertices$nepk,tolerance=1e-12)
    expect_equal(point$y,unname(GCEstim::NormEnt(best,model=FALSE)),tolerance=1e-12)
    if(kind=="TARWGCE")expect_equal(point$y,c(.03360347500108142,.99529587223169502,
      .44484511942118760,.94186252216217148,.99972801900330510),tolerance=1e-12)
    p <- entropy_plot(state)
    colours <- GCEstimjmv:::.gceTermColours(ids,jmvcore:::getGlobalTheme("default","jmv")$theme)
    expect_identical(p$built$data[[2]]$colour,unname(colours[point$predictor]))
    expect_identical(p$built$data[[1]]$colour,unname(colours[as.character(p$plot$data$predictor[order(p$plot$data$predictor,p$plot$data$support)])]))
    # A legacy plot4 omits a valid support; plot3 is the complete same-stage source.
    saved <- cache
    saved$plots$plot4$entropyComplete <- NULL
    saved$plots$plot4$data <- subset(saved$plots$plot4$data,support!=unique(state$data$support)[1])
    for(i in seq_along(saved$plots$plot4$overlays)) saved$plots$plot4$overlays[[i]]$predictor<-NULL
    before <- serialize(saved,NULL)
    local_mocked_bindings(cv.lmgce=function(...)stop("Unexpected CV"),
      lmgce=function(...)stop("Unexpected fit"),ridgetrace=function(...)stop("Unexpected ridge"),.package="GCEstim")
    for (pass in 1:2) {
      reopened <- entropy_analysis(kind)
      reopened$results$fitCache$.__enclos_env__$super$fromProtoBuf(state_message(unserialize(before)),character(),character())
      if(kind=="TARWGCE")
        reopened$results$plotRidge$.__enclos_env__$super$fromProtoBuf(state_message(a$results$plotRidge$state),character(),character())
      reopened$run()
      recovered <- GCEstimjmv:::.gcePlotImage(reopened$results$plot4,reopened$results$fitCache$state)$state
      expect_equal(nrow(recovered$data),n)
      expect_identical(recovered$data,cache$plots$plot3$data)
      expect_identical(serialize(reopened$results$fitCache$state,NULL),before)
      expect_identical(reopened$results$fitCache$state$coefficients,cache$coefficients)
      expect_identical(reopened$results$fitCache$state$ci,cache$ci)
      expect_identical(reopened$results$fitCache$state$fit,cache$fit)
      expect_identical(recovered$supportLines,state$supportLines)
      expect_identical(recovered$overlays[[1]]$y,point$y)
      expect_match(reopened$results[[".fitTime"]]$content,"Previous fit reused")
      expect_identical(entropy_plot(recovered)$built$data,p$built$data)
    }
  })
})

test_that("entropy legends wrap five and many long coefficient labels within the fixed width", {
  for(ids in list(c("(Intercept)","Z signal","A signal","F groupB","F groupC"),
                  paste("Long coefficient label with meaningful spaces and punctuation",1:10),
                  paste0("VeryLongUnbrokenCoefficientIdentifier",1:8))) {
    legend <- GCEstimjmv:::.gceEntropyLegend(ids)
    expect_identical(names(legend$labels),ids)
    expect_true(legend$columns>=1L)
    expect_lte(ceiling(length(ids)/legend$columns),5L)
    expect_true(all(nzchar(legend$labels)))
    d <- expand.grid(predictor=ids,support=c(1,2),stringsAsFactors=FALSE)
    d$nepk <- .7
    state <- list(which=4L,coefficientTerms=ids,data=d,overlays=list(),entropyComplete=TRUE)
    p <- entropy_plot(state)
    expect_identical(p$plot$guides$guides$colour$params$ncol,legend$columns)
    expect_true(p$plot$guides$guides$colour$params$theme$legend.byrow)
    expect_identical(p$plot$theme$legend.position,"bottom")
    expect_s3_class(p$plot$theme$legend.title,"element_blank")
    # Measure the actual legend on the fixed-size jamovi-compatible device.
    cap <- ragg::agg_capture(width=520,height=400,res=72)
    grob <- suppressMessages(ggplot2::ggplotGrob(p$plot))
    guide <- grob$grobs[[which(grob$layout$name=="guide-box-bottom")]]
    expect_lt(grid::convertWidth(sum(guide$widths),"in",valueOnly=TRUE),520/72)
    expect_lt(grid::convertHeight(sum(guide$heights),"in",valueOnly=TRUE),.65*400/72)
    dev.off()
  }
})

test_that("original and reopened validation states recover their actual missing rows", {
  # Captured read-only from the two named validation OMVs before this correction.
  states <- readRDS(test_path("_fixtures","entropy-legacy-states.rds"))
  local_mocked_bindings(cv.lmgce=function(...)stop("Unexpected CV"),
    lmgce=function(...)stop("Unexpected fit"),ridgetrace=function(...)stop("Unexpected ridge"),.package="GCEstim")
  for (origin in names(states)) for (kind in names(states[[origin]])) {
    cache <- states[[origin]][[kind]];before <- serialize(cache,NULL)
    expected <- if(kind=="TARWGCE")20L else 25L
    if(kind=="TARWGCE")expect_equal(nrow(cache$plots$plot4$data),15L)
    state <- GCEstimjmv:::.gcePlotImage(list(state=list(key="plot4",OLS=FALSE)),cache)$state
    expect_equal(nrow(state$data),expected)
    expect_identical(state$data,cache$plots$plot3$data)
    expect_identical(serialize(cache,NULL),before)
    expect_true(state$entropyComplete)
    p <- entropy_plot(state)
    expect_equal(nrow(p$built$data[[1]]),expected)
    point <- cache$plots$plot4$overlays[[1]]
    vertex <- merge(point,state$data,by.x=c("predictor","x"),by.y=c("predictor","support"))
    expect_equal(nrow(vertex),5L)
    expect_equal(vertex$y,vertex$nepk,tolerance=1e-12)
    expect_identical(p$built$data[[2]]$x,point$x)
    expect_identical(p$built$data[[2]]$y,point$y)
  }
  expect_identical(states$original,states$reopened)
})

test_that("an empty unrecoverable legacy entropy state renders an honest notice", {
  state <- GCEstimjmv:::.gceEntropyState(list(which=4L,data=data.frame(),overlays=list()))
  expect_equal(nrow(state$data),0L)
  p <- entropy_plot(state)
  expect_match(p$plot$labels$caption,"unavailable")
})

test_that("complete same-stage data can recover an entirely truncated entropy trajectory",{
  source <- list(which=3L,data=data.frame(support=c(1,2),predictor=c("Z","Z"),nepk=c(.2,.8)))
  state <- list(which=4L,coefficientTerms="Z",data=source$data[FALSE,],overlays=list())
  recovered <- GCEstimjmv:::.gceEntropyState(state,source)
  expect_identical(recovered$data,source$data)
  expect_true(recovered$entropyComplete)
  expect_null(recovered$entropyNotice)
})
