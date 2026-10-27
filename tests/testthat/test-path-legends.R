
# Exercise the same fixed-size export used by the existing-analysis callbacks.
path_legend_plot <- function(renderer, state) {
  fun <- get(renderer,asNamespace("GCEstimjmv"))
  env <- new.env(parent=environment(fun))
  env$print <- function(x){env$p<-x;invisible(x)};environment(fun)<-env
  g<-jmvcore:::getGlobalTheme("default","jmv")
  expect_true(suppressMessages(fun(list(state=state),g$ggtheme,g$theme)))
  env$p
}

test_that("coefficient-path legends wrap without changing geometry, colours or extrema",{
  cases <- list(c("(Intercept)","Z signal","\u0060F group\u0060B","F groupB","A signal"),
    paste("Long coefficient name with meaningful punctuation (factor level)",1:10),
    paste0("VeryLongUnbrokenCoefficientIdentifier",1:8))
  for(ids in cases)for(renderer in c(".gceRenderRidge",".gceRenderDiagnostic")) {
    ids <- rev(ids)
    d <- expand.grid(id=ids,x=c(1,2,3),KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
    d$y<-rep(seq_along(ids),3)+rep(c(-.1,0,.1),each=length(ids))
    if(renderer==".gceRenderRidge") {
      state<-list(coefficientTerms=ids,data=data.frame(variable=d$id,lambda=d$x,coefficient=d$y),
        extrema=data.frame(variable=ids,max_x=3,max_y=seq_along(ids)+.1,min_x=1,min_y=seq_along(ids)-.1))
    } else state<-list(which=3L,coefficientTerms=ids,
      data=data.frame(predictor=d$id,support=d$x,estimate=d$y),overlays=list(),
      supportLines=list(min=1,oneSE=2,manual=3))
    before<-serialize(state,NULL)
    p<-path_legend_plot(renderer,unserialize(before));built<-ggplot2::ggplot_build(p)
    expect_identical(serialize(state,NULL),before)
    expect_identical(p$theme$legend.position,"bottom")
    expect_s3_class(p$theme$legend.title,"element_blank")
    legend<-GCEstimjmv:::.gceCoefficientLegend(ids)
    expect_identical(p$guides$guides$colour$params$ncol,legend$columns)
    expect_true(p$guides$guides$colour$params$theme$legend.byrow)
    expect_identical(unname(built$plot$scales$get_scales("colour")$get_labels(ids)),unname(legend$labels[ids]))
    colours<-GCEstimjmv:::.gceTermColours(ids,jmvcore:::getGlobalTheme("default","jmv")$theme)
    expect_identical(unname(built$plot$scales$get_scales("colour")$map(ids)),unname(colours))
    expect_equal(sort(built$data[[1]]$y),sort(d$y))
    expect_equal(sort(built$data[[1]]$x),sort(d$x))
    if(renderer==".gceRenderRidge") {
      expect_identical(built$data[[2]]$shape,rep(24,length(ids)))
      expect_identical(built$data[[3]]$shape,rep(25,length(ids)))
      expect_identical(built$data[[2]]$size,rep(1,length(ids)))
      expect_identical(built$data[[3]]$size,rep(1,length(ids)))
      expect_equal(built$data[[2]]$x,state$extrema$max_x)
      expect_equal(built$data[[2]]$y,state$extrema$max_y)
      expect_equal(built$data[[3]]$x,state$extrema$min_x)
      expect_equal(built$data[[3]]$y,state$extrema$min_y)
      expect_identical(built$data[[2]]$colour,unname(colours[ids]))
    } else {
      refs<-Filter(function(x)"xintercept"%in%names(x),built$data)
      expect_identical(vapply(refs,function(x)x$xintercept[1],numeric(1)),c(1,2,3))
      expect_identical(vapply(refs,function(x)x$linetype[1],character(1)),c("dotted","dashed","dotted"))
      expect_identical(vapply(refs,function(x)x$colour[1],character(1)),c("black","black","green4"))
    }
    cap<-ragg::agg_capture(width=520,height=400,res=72)
    grob<-suppressMessages(ggplot2::ggplotGrob(p))
    guide<-grob$grobs[[which(grob$layout$name=="guide-box-bottom")]]
    expect_lt(grid::convertWidth(sum(guide$widths),"in",valueOnly=TRUE),520/72)
    expect_lt(grid::convertHeight(sum(guide$heights),"in",valueOnly=TRUE),.7*400/72)
    dev.off()
    file<-tempfile(fileext=".png")
    suppressMessages(ggplot2::ggsave(file,p,device=ragg::agg_png,width=520/72,height=400/72,dpi=72,bg="white"))
    expect_gt(file.info(file)$size,0)
    expect_identical(dim(png::readPNG(file))[1:2],c(400L,520L))
    unlink(file)
  }
})
