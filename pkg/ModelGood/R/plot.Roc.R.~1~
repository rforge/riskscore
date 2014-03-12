plot.Roc <- function(x,
                     ylab="True detection rate",
                     xlab="False positive rate",
                     models,
                     type="s",
                     roc=TRUE,
                     shadow=FALSE,
                     simu=FALSE,
                     BCV=FALSE,
                     apparent=FALSE,
                     noinf=FALSE,
                     control,
                     grid=FALSE,
                     diag=FALSE,
                     box=TRUE,
                     lwd=2,
                     lty,
                     col,
                     add=FALSE,
                     axes=TRUE,
                     las=2,
                     legend=TRUE,
                     percent=TRUE,
                     ...
                     ){
  
  # find the models 
  # --------------------------------------------------------------------
  if (missing(models)){
    found <-  1:length(x$models)
    models <- names(x$models)
  }
  else{
    if (!is.numeric(models)){
      stopifnot(all(found <- match(models,names(x$models),nomatch=0))>0)
      models <- names(x$models)[found]
    }
    else{
      stopifnot(all(found <- match(models,1:length(x$models),nomatch=0))>0)
      models <- names(x$models)[found]
    }
  }

  # find the lines to be plotted 
  # --------------------------------------------------------------------
  
  do <- list("Roc"=0,"AppRoc"=0,"BootcvRoc"=0,"NoInfRoc"=0)
  if (roc==TRUE) do[["Roc"]] <- 1 
  if (apparent==TRUE) do[["AppRoc"]] <- do[["Roc"]] + 1
  if (BCV==TRUE) do[["BootcvRoc"]] <- do[["Roc"]] + (do[["AppRoc"]]>0) + 1
  if (noinf==TRUE) do[["NoInfRoc"]] <- do[["Roc"]] + (do[["AppRoc"]]>0) + (do[["BootcvRoc"]]>0) + 1
  
  # set default col and lty
  # --------------------------------------------------------------------
  
  n.models <- length(x$models)
  def.col <- c("Roc"=1,"AppRoc"=1,"BootcvRoc"=1,"NoInfRoc"=1)
  def.lty <- rep(1,n.models)
  names(def.lty) <- models
  
  if (missing(col)){
    col <- as.list(1:n.models)
  }
  else{
    tmpCol <- as.list(1:n.models)
    whoM <- match(models,names(x$models),nomatch=FALSE)
    tmpCol[whoM] <- col[1:length(whoM>0)]
    col <- tmpCol
  }
  names(col) <- names(x$models)
  if (missing(col))
  col <- lapply(col,function(x){def.col*x})
  
  if (missing(lty)){
    #    "Roc" 1
    #    "AppRoc" 2
    #    "BCVRoc" 3
    #    "NoInf" 4
    lty <- lapply(col,function(x){
      cumsum(unlist(do))
    })
  }
  
  if (missing(lwd)) lwd <- 2
  
  legend.DefaultArgs <- list(legend=names(x$models),
                             lwd=sapply(lwd,function(x)x[1]),
                             col=sapply(col,function(x)x[1]),
                             lty=sapply(lty,function(x)x[1]),
                             cex=1.5,
                             bty="n",
                             y.intersp=1.3,
                             x="bottomright")
  
  plot.DefaultArgs <- list(x=0,
                           y=0,
                           type = "n",
                           ylim = 0:1,
                           xlim = 0:1,
                           xlab = xlab,
                           ylab = ylab)
  
  xaxis.DefaultArgs <- list(at=seq(0,1,.25),pos=0)
  yaxis.DefaultArgs <- list(at=seq(0,1,.25),pos=0)
  smartA <- MgResolveSmartArgs(call=  list(...),
                               keys=c("plot","legend","xaxis","yaxis"),
                               ignore=c("x","add","col","lty","lwd","ylim","xlim","xlab","ylab","legend","percent","grid","box","axes","roc","type","shadow","simu","BCV","apparent","noinf","control","diag"),
                               defaults=list("plot"=plot.DefaultArgs,
                                 "legend"=legend.DefaultArgs,
                                 "xaxis"=xaxis.DefaultArgs,
                                 "yaxis"=yaxis.DefaultArgs),
                               forced=list("plot"=list(axes=FALSE),
                                 "xaxis"=list(side=1),
                                 "yaxis"=list(side=2)),
                               verbose=TRUE)
  
  # plot an empty frame
  # --------------------------------------------------------------------
  if (!add) {
    do.call("plot",smartA$plot)
    if (axes){
      if (percent & is.null(smartA$xaxis$labels))
        smartA$xaxis$labels <- paste(100*smartA$xaxis$at,"%")
      do.call("axis",smartA$xaxis)
      if (las==2){
        oldlas=par()$las
        par(las=2)}
      if (percent & is.null(smartA$yaxis$labels))
        smartA$yaxis$labels <- paste(100*smartA$yaxis$at,"%")
      do.call("axis",smartA$yaxis)
      if (las==2){par(las=oldlas)}
    }
    if (grid==TRUE) abline(h = 0:10/10, v = 0:10/10, col = gray(0.9))
    if (diag==TRUE) abline(0, 1, col = gray(0.4))
    if (box==TRUE)  box()
  }
  
  # adding the lines
  # --------------------------------------------------------------------
  ## produce a shadow for the roc curve of each model
  ## showing the performance in the bootstrap runs
  
  if (!(is.logical(shadow) && shadow==FALSE)){
    ## ?hcl different colors for shading
    if (is.logical(shadow)) shadow.control <- list()
    else shadow.control <- shadow
    nix <- lapply(found,function(m){
      nix <- lapply(x$BootcvRocList[[m]],function(x){
        shadow.control <- c(shadow.control,list(type=type,col="gray77"))
        shadow.control <- shadow.control[!duplicated(names(shadow.control))]
        do.call("lines.default",c(list(1-x$Specificity, x$Sensitivity),shadow.control))
      })
    })
  }
  if (!(is.logical(simu) && simu==FALSE)){
    ## ?hcl different colors for shading
    if (is.logical(simu)) simu.control <- list()
    else simu.control <- simu
    nix <- lapply(found,function(m){
      nix <- lapply(x$NoInfRocList[[m]],function(x){
        simu.control <- c(simu.control,list(type=type,col="gray77"))
        simu.control <- simu.control[!duplicated(names(simu.control))]
        do.call("lines.default",c(list(1-x$Specificity, x$Sensitivity),simu.control))
      })
    })
  }

  
  if (missing(control))
    control <- do
  
  for (w in 1:4){
    W <- names(do)[w]
    if (do[[w]]>0){
      for (m in found){
        w.control <- control[[w]]
        w.control <- c(w.control,list(type=type,lwd=lwd,lty=lty[[m]][w],col=col[[m]][w]))
        w.control <- w.control[!duplicated(names(w.control))]
        X <- 1-x[[W]][[m]]$Specificity
        Y <- x[[W]][[m]]$Sensitivity
        do.call("lines.default",c(list(X,Y), w.control))
      }
    }
  }
  # legend
  # --------------------------------------------------------------------
  if(legend==TRUE && !add && !is.null(names(x$models))){
    save.xpd <- par()$xpd
    par(xpd=TRUE)
    do.call("legend",smartA$legend)
    par(xpd=save.xpd)
  }
}

lines.Roc <- function(x,...){
  plot.Roc(x,add=TRUE,...)
}
