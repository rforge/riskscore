calibrationPlot <- function(object,
                            formula,
                            data,
                            method="snne",
                            bandwidth=NULL,
                            smooth=TRUE,
                            verbose=FALSE,
                            add=FALSE,
                            showData=(!add),
                            diag=TRUE,
                            legend=TRUE,
                            axes=TRUE,
                            background=FALSE,
                            Grid=background,
                            xlim,
                            ylim,
                            xlab,
                            ylab,
                            col,
                            lwd,
                            lty,
                            pch,
                            cause=1,
                            percent=TRUE,
                            ...){
  
  # {{{ find number of objects and lines
  cobj=class(object)
  if (cobj[[1]]!="list"){
    object <- list(object)
  }
  if (is.null(names(object))){
    names(object) <- sapply(object,function(o)class(o)[1])
  }
  else{
    names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
  }
  names(object) <- make.names(names(object),unique=TRUE)
  nlines <- length(object)
  # }}}
  # {{{ lines types
  if (missing(lwd)) lwd <- rep(3,nlines)
  if (missing(col)) col <- 1:nlines
  if (missing(lty)) lty <- rep(1, nlines)
  if (missing(pch)) pch <- rep(1, nlines)
  if (length(lwd) < nlines) lwd <- rep(lwd, nlines)
  if (length(lty) < nlines) lty <- rep(lty, nlines)
  if (length(col) < nlines) col <- rep(col, nlines)
  if (length(pch) < nlines) pch <- rep(pch, nlines)
  # }}}
  # {{{ data & formula
  if (missing(data)){
    data <- eval(object[[1]]$call$data)
    if (match("data.frame",class(data),nomatch=0)==0)
      stop("Argument data is missing.")
    else
      if (verbose)
        warning("Argument data is missing. I use the data from the call to the first model instead.")
  }
  if (missing(formula)){
    formula <- eval(object[[1]]$call$formula)
    if (match("formula",class(formula),nomatch=0)==0)
      stop("Argument formula is missing.")
    else if (verbose)
      warning("Argument formula is missing. I use the formula from the call to the first model instead.")
  }
  m <- model.frame(formula,data,na.action=na.fail)
  Y <- model.response(m)
  # }}}
  # {{{ call smartControls
  axis1.DefaultArgs <- list(side=1,las=1)
  axis2.DefaultArgs <- list(side=2,las=2)
  legend.DefaultArgs <- list(legend=names(object),
                             lwd=lwd,
                             col=col,
                             lty=lty,
                             cex=1.5,
                             bty="n",
                             y.intersp=1.3,
                             x="topleft")
  Grid.DefaultArgs <- list(vertical=NULL,horizontal=NULL,col=c("gray88","gray99"))
  background.DefaultArgs <- list(col="gray77")
  lines.DefaultArgs <- list(type="l")
  if (missing(ylim)){
    ylim <- c(0,1)
  }
  if (missing(xlim)){
    xlim <- c(0,1)
  }
  plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,xlab = "Predicted probabilities",ylab = "Observed probabilities")
  smartA <- prodlim:::SmartControl(call= list(...),
                                   keys=c("plot","lines","legend","Grid","background","axis1","axis2"),
                                   ignore=NULL,
                                   ignore.case=TRUE,
                                   defaults=list("plot"=plot.DefaultArgs,"lines"=lines.DefaultArgs,"legend"=legend.DefaultArgs,"Grid"=Grid.DefaultArgs,"background"=background.DefaultArgs,"axis1"=axis1.DefaultArgs,"axis2"=axis2.DefaultArgs),
                                   forced=list("plot"=list(axes=FALSE),"axis1"=list(side=1)),
                                   verbose=TRUE)
  # }}}
  # {{{ plot an empty frame
  if (add==FALSE){
    do.call("plot",smartA$plot)
    if (axes){
      if (percent){
        if (is.null(smartA$axis2$at)){
          smartA$axis2$at <- seq(0,1,.25)
          smartA$axis2$labels <- paste(100*smartA$axis2$at,"%")
        }
        else{
          smartA$axis2$labels <- paste(100*smartA$axis2$at,"%")
        }
        if (is.null(smartA$axis1$at)){
          smartA$axis1$at <- seq(0,1,.25)
          smartA$axis1$labels <- paste(100*smartA$axis1$at,"%")
        }
        else{
          smartA$axis1$labels <- paste(100*smartA$axis1$at,"%")
        }
      }
      do.call("axis",smartA$axis2)
      do.call("axis",smartA$axis1)
      do.call("axis",smartA$axis2)
    }
    if (background)
      do.call("prodlim:::background",smartA$background)
    if (Grid){
      if (is.null(smartA$Grid$horizontal) && !is.null(smartA$axis2$at))
        smartA$Grid$horizontal <- smartA$axis2$at
      if (is.null(smartA$Grid$vertical) && !is.null(smartA$axis1$at))
        smartA$Grid$vertical <- smartA$axis1$at
      do.call("prodlim:::Grid",smartA$Grid)
    }
  }
  if (diag & add==FALSE)
    segments(x0=0,y0=0,x1=1,y1=1,col="gray77",lwd=2,xpd=FALSE)
  ##   do.call("abline",c(list(a=0,b=1),list(col="gray77",lwd=2,xpd=FALSE)))
  
  # }}}
  # {{{ add one smoothed calibration line for each model
  method <- match.arg(method,c("fixed","snne","loess","kernel"))
  predictions <- lapply(1:nlines,function(f){
    p=do.call("predictStatusProb",list(object[[f]],newdata=data))
    ## if (class(object[[f]])=="matrix") p <- p[neworder,]
    if (showData) points(p,Y,col="gray",xpd=NA)
    switch(method,
           "fixed"={
             x <- unique(c(quantile(p,seq(0,1,0.1))))
             x[1] <- 0
             ybar <- sapply(2:length(x),function(q){
               mean(Y[p>x[q-1] & p<=x[q]])
             })
             mid <- x[-length(x)]+diff(x)/2
             if (smooth)
               lines(lowess(cbind(mid,ybar)),col=col[f],lty=lty[f],lwd=lwd[f])
             else 
               lines(mid,ybar,col=col[f],lty=lty[f],lwd=lwd[f])
             
             ## y <- sapply(unique(p),function(prob){
             ## nbh <- (1:NROW(p)) prob
             ## jack
             ## }
             ## browser()
           },
           "snne"={
             plotFrame=prodlim:::meanNeighbors(x=p,y=as.numeric(as.character(Y)),bandwidth=bandwidth)
             lines(plotFrame$uniqueX,
                   plotFrame$averageY,
                   col=col[f],
                   lty=lty[f],
                   lwd=lwd[f])
           },
           "loess"={
             loess.args=NULL
             loess.DefaultArgs <- list(family="symmetric",control=loess.control(iterations=0),span=.3)
             loess.args <- c(formula=Y~pp,loess.args,loess.DefaultArgs)
             loess.args <- loess.args[!duplicated(names(loess.args))]
             pp=sort(p)
             smoothY=do.call("loess",loess.args)
             plotFrame=data.frame(x=pp,y=smoothY$fitted)
             with(na.omit(plotFrame),lines(x,y,col=col[f],lwd=lwd[f],lty=lty[f]))
           },
           "kernel"={
             ##              require(kknn)
             message("kernel smoothing not yet implemented")
           })
    # lines(pseudoFrame$uniqueX,smooth(pseudoFrame$averageY),col=col[f],lwd=lwd)
    c(p)
  })
  # }}}
  # {{{ legend
  ## if (missing(legend)) legend=ifelse(length(object)==1,FALSE,TRUE)
  ## if (missing(legend.legend)) legend.legend=names(object)
  if (add==FALSE & legend==TRUE){
    do.call("legend",smartA$legend)
  }
  ## if (legend)
  ## legend(0,1,legend=legend.legend,lwd=lwd,col=col,bty="n")
  # }}}
  # {{{ invisibly output the response and the predictions
  out <- list(Y=Y,predictions=predictions)
  invisible(out)
  # }}}
}
