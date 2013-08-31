calPlot <- function(object,
                    predTime,
                    formula,
                    data,
                    splitMethod="none",
                    B=1,
                    M,
                    giveToModel=NULL,
                    method="nne",
                    q=10,
                    outcome=c("pseudo","prodlim"),
                    bandwidth=NULL,
                    verbose=FALSE,
                    add=FALSE,
                    showPseudo=ifelse(add||(outcome!="pseudo"),FALSE,TRUE),
                    jack.density=55,
                    diag=ifelse(add,FALSE,TRUE),
                    legend=ifelse(add,FALSE,TRUE),
                    axes=ifelse(add,FALSE,TRUE),
                    background=FALSE,
                    Grid=background,
                    xlim,
                    ylim,
                    xlab = "Predicted probabilities",
                    ylab = "Observed probabilities",
                    col,
                    lwd,
                    lty,
                    pch,
                    cause=1,
                    percent=TRUE,
                    na.action=na.fail,
                    cores=1,
                    ...){
  
  # {{{ find number of objects and lines
  cobj=class(object)[[1]]
  if (cobj!="list"){
    object <- list(object)
  }
  if (is.null(names(object))){
    names(object) <- sapply(object,function(o)class(o)[1])
  }
  else{
    names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
  }
  names(object) <- make.names(names(object),unique=TRUE)
  NF <- length(object)

  # }}}
  # {{{ lines types
  if (missing(lwd)) lwd <- rep(3,NF)
  if (missing(col)) col <- 1:NF
  if (missing(lty)) lty <- rep(1, NF)
  if (missing(pch)) pch <- rep(1, NF)
  if (length(lwd) < NF) lwd <- rep(lwd, NF)
  if (length(lty) < NF) lty <- rep(lty, NF)
  if (length(col) < NF) col <- rep(col, NF)
  if (length(pch) < NF) pch <- rep(pch, NF)
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
  m <- model.frame(formula,data,na.action=na.action)
  response <- model.response(m)
  if (match("Surv",class(response),nomatch=FALSE))
    model.type <- "survival"
  else
    model.type <- attr(response,"model")
  neworder <- order(response[,"time"],-response[,"status"])
  response <- response[neworder,,drop=FALSE]
  Y <- response[,"time"]
  status <- response[,"status"]
  data <- data[neworder,]
  # }}}
  # {{{ prediction timepoint 
  
  if (missing(predTime))
    predTime <- median(Y)
  else
    if (length(predTime)>1)
      stop("Please specify only one time point predTime.")

  # }}}
  # {{{ compute pseudo-values
  #  require(pseudo)
  #  jack=pseudosurv(time=Y,event=status,tmax=predTime)[[3]]
  if (model.type=="competing.risks"){
    predictHandlerFun <- "predictEventProb"
  }
  else{
    predictHandlerFun <- "predictSurvProb"
  }
  margForm <- reformulate("1",response=formula[[2]])
  margFit <- prodlim(margForm,data=data)
  jack <- jackknife(margFit,cause=cause,times=predTime)

  # }}}
  # {{{ call smartControls

  axis1.DefaultArgs <- list(side=1,las=1)
  axis2.DefaultArgs <- list(side=2,las=2)
  legend.DefaultArgs <- list(legend=names(object),lwd=lwd,col=col,lty=lty,cex=1.5,bty="n",y.intersp=1.3,x="topleft")
  Grid.DefaultArgs <- list(vertical=NULL,horizontal=NULL,col=c("gray88","gray99"))
  background.DefaultArgs <- list(col="gray77")
  lines.DefaultArgs <- list(type="l")
  if (missing(ylim)){
    if (showPseudo)
      ylim <- c(min(jack),max(jack))
    else
      ylim <- c(0,1)
  }
  if (missing(xlim)){
    xlim <- c(0,1)
  }
  plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,ylab=ylab,xlab=xlab)
  smartA <- prodlim:::SmartControl(call= list(...),keys=c("plot","lines","legend","Grid","background","axis1","axis2"),ignore=NULL,ignore.case=TRUE,defaults=list("plot"=plot.DefaultArgs,"lines"=lines.DefaultArgs,"legend"=legend.DefaultArgs,"Grid"=Grid.DefaultArgs,"background"=background.DefaultArgs,"axis1"=axis1.DefaultArgs,"axis2"=axis2.DefaultArgs),forced=list("plot"=list(axes=FALSE),"axis1"=list(side=1)),verbose=TRUE)

  # }}}
  # {{{ splitmethod
  splitMethod <- resolvesplitMethod(splitMethod=splitMethod,B=B,N=NROW(data),M=M)
  k <- splitMethod$k
  B <- splitMethod$B
  N <- splitMethod$N
  NF <- length(object) 
  # }}}
  # {{{ cv, predictions and expectations
  # {{{ ---------------------------Apparent predictions---------------------------
  apppred <- do.call("cbind",lapply(1:NF,function(f){
    if (class(object[[f]][[1]])=="matrix"){
      apppred <- object[[f]][[1]][neworder,]
    } else{
      if (model.type=="competing.risks"){
        apppred <- do.call(predictHandlerFun,list(object[[f]],newdata=data,times=predTime,cause=cause))
      }
      else{
        apppred <- do.call(predictHandlerFun,list(object[[f]],newdata=data,times=predTime))
      }
    }
  }))
  colnames(apppred) <- names(object)
  apppred <- data.frame(jack=jack,apppred)
  if (splitMethod$internal.name %in% c("noPlan")){
    predframe <- apppred
  }
  # }}}
  # {{{--------------k-fold and leave-one-out CrossValidation-----------------------
  if (splitMethod$internal.name %in% c("crossval","loocv")){
    groups <- splitMethod$index[,1,drop=TRUE]
    cv.list <- lapply(1:k,function(g){
      if (verbose==TRUE) internalTalk(g,k)
      id <- groups==g
      train.k <- data[!id,,drop=FALSE]
      val.k <- data[id,,drop=FALSE]
      model.pred <- lapply(1:NF,function(f){
        extraArgs <- giveToModel[[f]]
        fit.k <- internalReevalFit(object=object[[f]],data=train.k,step=paste("CV group=",k),silent=FALSE,verbose=verbose)
        if (predictHandlerFun == "predictEventProb"){      
          do.call(predictHandlerFun,c(list(object=fit.k,newdata=val.k,times=predTime,cause=cause,train.data=train.k),extraArgs))
        }
        else{
          do.call(predictHandlerFun,c(list(object=fit.k,newdata=val.k,times=predTime,train.data=train.k),extraArgs))
        }
      })
      model.pred
    })
    predframe <- do.call("cbind",lapply(1:NF,function(f){
      pred <- do.call("rbind",lapply(cv.list,function(x)x[[f]]))
      if (splitMethod$internal.name!="loocv"){
        pred <- pred[order(order(groups)),]
      }
      pred
    }))
    colnames(predframe) <- names(object)
    predframe <- cbind(data.frame(jack=jack),predframe)
  }
  # }}} 
  # {{{ ----------------------BootstrapCrossValidation----------------------
  
  if (splitMethod$internal.name %in% c("Boot632plus","BootCv","Boot632")){
    if (splitMethod$internal.name %in% c("Boot632plus","Boot632")){
      stop("Don't know how to do the 632(+) for the calibration curve.")
    }
    ResampleIndex <- splitMethod$index
    ## predframe <- do.call("rbind",lapply(1:B,function(b){
    ## predframe <- matrix
    pred.list <- mclapply(1:B,function(b){
      if (verbose==TRUE) internalTalk(b,B)
      jackRefit <- FALSE
      vindex.b <- match(1:N,unique(ResampleIndex[,b]),nomatch=0)==0
      val.b <- data[vindex.b,,drop=FALSE]
      if (jackRefit){
        margFit.b <- prodlim(margForm,data=val.b)
        jack.b <- jackknife(margFit.b,cause=cause,times=predTime)
      }
      else{
        jack.b <- jack[match(1:N,unique(ResampleIndex[,b]),nomatch=0)==0]
      }
      train.b <- data[ResampleIndex[,b],,drop=FALSE]
      frame.b <- data.frame(jack=jack.b)
      bootpred <- do.call("cbind",lapply(1:NF,function(f){
        fit.b <- internalReevalFit(object=object[[f]],data=train.b,step=b,silent=FALSE,verbose=verbose)
        extraArgs <- giveToModel[[f]]
        if (predictHandlerFun == "predictEventProb"){
          try2predict <- try(pred.b <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=val.b,times=predTime,cause=cause,train.data=train.b),extraArgs)))
        }
        else{
          try2predict <- try(pred.b <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=val.b,times=predTime,train.data=train.b),extraArgs)))
        }
        rm(fit.b)
        gc()
        if (inherits(try2predict,"try-error")==TRUE){
          rep(NA,NROW(val.b))
        }
        pred.b
      }))
      colnames(bootpred) <- names(object)
      cbind(frame.b,bootpred)
    },mc.cores=cores)
    predframe <- do.call("rbind",pred.list)
    rm(pred.list)
  }

  # }}}
  # }}}
  # {{{ smoothing
  method <- match.arg(method,c("quantile","nne","loess"))
  outcome <- match.arg(outcome,c("pseudo","prodlim"))
  plotFrames <- lapply(1:NF,function(f){
      p <- predframe[,f+1]
      jackF <- predframe[,1]
      switch(method,
             "quantile"={
                 groups <- quantile(p,seq(0,1,1/q))
                 xgroups <- (groups[-(q+1)]+groups[-1])/2
                 if (outcome=="pseudo"){
                     plotFrame=data.frame(x=xgroups,y=tapply(jackF,cut(p,groups,include.lowest=TRUE),mean))
                 }
                 else{
                     pcut <- cut(p,groups,include.lowest=TRUE)
                     form.pcut <- reformulate("pcut",response=formula[[2]])
                     y <- unlist(predict(prodlim(form.pcut,data=cbind(data,p=pcut)),
                                         cause=cause,
                                         newdata=data.frame(pcut=levels(pcut)),
                                         times=predTime,
                                         type=ifelse(model.type=="competing.risks","cuminc","surv")))
                     plotFrame=data.frame(x=xgroups,y=y)
                 }
             },
             "nne"={
                 if (outcome=="pseudo"){
                     ## Round probabilities to 2 digits
                     ## to avoit memory explosion ...
                     ## a difference in the 3 digit should
                     ## not play a role for the patient.
                     p <- round(p,2)
                     bw <- neighborhood(apppred[,f+1])$bandwidth
                     ## print(bw)
                     nbh <- meanNeighbors(x=p,y=jackF,bandwidth=bw)
                     plotFrame <- data.frame(x=nbh$uniqueX,y=nbh$averageY)
                 }else{
                     form.p <- reformulate("p",response=formula[[2]])
                     y <- unlist(predict(prodlim(form.p,data=cbind(data,p=p)),
                                         cause=cause,
                                         newdata=data.frame(p=sort(p)),
                                         times=predTime,
                                         type=ifelse(model.type=="competing.risks","cuminc","surv")))
                     plotFrame <- data.frame(x=sort(p),y=y)
                 }
             },
             "loess"={
                 loess.args=NULL
                 loess.DefaultArgs <- list(family="symmetric",control=loess.control(iterations=0),span=.3)
                 pp <- sort(p)
                 lframe <- data.frame(pp=pp,jack=jackF)
                 loess.args <- c(formula=jack~pp,data=lframe,loess.args,loess.DefaultArgs)
                 loess.args <- loess.args[!duplicated(names(loess.args))]
                 smoothJack <- do.call("loess",loess.args)
                 plotFrame <- data.frame(x=pp,y=smoothJack$fitted)
             })
  })
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
  if (diag){
      segments(x0=0,y0=0,x1=1,y1=1,col="gray77",lwd=2,xpd=FALSE)
  }
  ##   do.call("abline",c(list(a=0,b=1),list(col="gray77",lwd=2,xpd=FALSE)))
  # }}}
  # {{{ add lines and pseudovalues
  nix <- lapply(1:NF,function(f){
      plotFrame <- plotFrames[[f]]
      with(na.omit(plotFrame),lines(x,y,col=col[f],lwd=lwd[f],lty=lty[f],type=ifelse(method=="quantile","b","l")))
      ccrgb=as.list(col2rgb(col[f],alpha=TRUE))
      names(ccrgb) <- c("red","green","blue","alpha")
      ccrgb$alpha <- jack.density
      jack.col <- do.call("rgb",c(ccrgb,list(max=255)))
      if (showPseudo) {
          points(apppred[,f+1],apppred[,1],col=jack.col)
      }
  })

  # }}}
  # {{{ legend
  ## if (missing(legend)) legend=ifelse(length(object)==1,FALSE,TRUE)
  ## if (missing(legend.legend)) legend.legend=names(object)
  if(legend){
    do.call("legend",smartA$legend)
  }
  ## if (legend)
  ## legend(0,1,legend=legend.legend,lwd=lwd,col=col,bty="n")

  # }}}
  # {{{ invisibly output the jackknife pseudo-values
  out <- list(predTime=predTime,pseudoFrame=predframe)
  invisible(out)
  # }}}
}
