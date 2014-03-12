
# {{{ UseMethod

Roc <- function(object,...){
  UseMethod("Roc",object=object)
}

# }}}
# {{{ Roc.default, Roc.glm,etc.

Roc.default <- function(object,
                        y,
                        breaks,
                        crRatio=1,
                        ## pv=FALSE,
                        ## confint=FALSE,
                        ## confint.method="exact",
                        keep.tables=FALSE,
                        keep.breaks=FALSE){
  N <- length(object)
  if(length(y)!=N) stop("Arguments must have the same length")
  if(length(unique(y))!=2) stop("y must be binary")
  Disease <- as.integer(as.character(factor(y,labels=c("0","1"))))
  count.DiseasePos <- sum(Disease==1)
  count.DiseaseNeg <- sum(Disease==0)
  ## print(breaks)
  if (missing(breaks) || is.null(breaks))
    breaks <- sort(unique(object))
  else
    breaks <- sort(unique(breaks))
  if (length(breaks)>1 & !is.factor(object)){
    eval.times <- breaks- min(diff(breaks))/2
    count.TestPos <- N-MgSindex2(jump.times=object,eval.times=eval.times)
    count.TestNeg <- N-count.TestPos
    a <- count.DiseasePos-MgSindex2(jump.times=object[Disease==1],eval.times=eval.times)
    b <- count.DiseaseNeg-MgSindex2(jump.times=object[Disease==0],eval.times=eval.times)
  }
  else{
    tabx <- table(object)
    count.TestPos <- tabx
    count.TestNeg <- tabx
    tabxpos <- table(object[Disease==1])
    tabxneg <- table(object[Disease==0])
    a <- count.DiseasePos-tabxpos
    b <- count.DiseaseNeg-tabxneg
  }
  c <- count.DiseasePos-a
  d <- count.DiseaseNeg-b
  sens <- crRatio*a/count.DiseasePos
  spec <- d/count.DiseaseNeg
  out <- list("Sensitivity"=c(sens,0),"Specificity"=c(spec,1))
  ##   if (confint==TRUE){
  ##     tmp.sens <- binconf(x=a,n=count.DiseasePos,method=confint.method)
  ##     sens <- tmp.sens[,"PointEst"]
  ##     ci.sens <- tmp.sens[,c("Lower","Upper")]
  ##     tmp.spec <- binconf(x=d,n=count.DiseaseNeg,method=confint.method)
  ##     spec <- tmp.spec[,"PointEst"]
  ##     ci.spec <- tmp.spec[,c("Lower","Upper")]
  ##     out <- list("Sensitivity"=c(sens,0),"Specificity"=c(spec,1),
  ##                 "CI.Sens"=rbind(ci.sens,c(0,0)),"CI.Spec"=rbind(ci.spec,c(1,1)))
  ##   }
  ##   if (pv==TRUE){
  ##     if (confint==TRUE){
  ##       tmp.ppv <- binconf(x=a,n=count.TestPos,method=confint.method)
  ##       ppv <- tmp.ppv[,"PointEst"]
  ##       ci.ppv <- tmp.ppv[,c("Lower","Upper")]
  ##       tmp.npv <- binconf(x=d,n=count.TestNeg,method=confint.method)
  ##       npv <- tmp.npv[,"PointEst"]
  ##       ci.npv <- tmp.npv[,c("Lower","Upper")]
  ##       out <- c(out,list("PPV"=c(ppv,1),"NPV"=c(npv,npv[length(npv)]),
  ##                         list("CI.PPV"=rbind(ci.ppv,c(1,1)),
  ##                              "CI.NPV"=rbind(ci.npv,ci.npv[NCOL(ci.npv),]))))
  ##     }
  ##     else{
  ##   ppv <- a/count.TestPos
  ##   npv <- d/count.TestNeg
  ## }
  ##   out <- c(out,list("PPV"=c(ppv,1),"NPV"=c(npv,npv[length(npv)])))
  ## }
  if (keep.breaks==TRUE)
    out <- c(out,list(breaks=breaks))
  if (keep.tables==TRUE){
    out <- c(out,list(tables=data.frame(a,b,c,d)))
  }
  ##   if (confint==TRUE)
  ##     out$confint.method <- confint.method
  ## class(out) <- "Roc"
  out
}


Roc.formula <- function(object,formula,data,crRatio=1,...){
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m <- m[match(c("","formula","data","subset","na.action"), names(m), nomatch = 0)]
  if (missing(data)) Terms <- terms(formula)
  else Terms <- terms(formula,data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m1 <- eval(m, parent.frame())
  if (NROW(m1) == 0) stop("No (non-missing) status values")
  m$formula <- delete.response(Terms)
  m2 <- eval(m, parent.frame())
  if (NROW(m2) == 0) stop("No (non-missing) test values")
  y <- model.extract(m1,"response")
  x <- m2[,1,drop=TRUE]
  if(any(match(c("Surv","Hist"),all.names(formula),nomatch=0)>0)){
    stop("Survival models not supported. Refer to the R package `pec'.")
  }
  else{
    out <- list("Roc"=list(threshModel=Roc.default(object,y, crRatio=crRatio,...)),
                call=match.call(),
                models=c("threshModel"=formula))
    class(out) <- "Roc"
    out
  }
}

Roc.glm <- function(object,formula,data,...){
  stopifnot(object$family$family=="binomial")
  Roc.list(object=list(object),formula,data,...)
}

Roc.lrm <- function(object,formula,data,...){
  Roc.list(object=list(object),formula,data,...)
}


Roc.rpart <- function(object,formula,data,...){
  Roc.list(object=list(object),formula,data,...)
}

Roc.randomForest <- function(object,formula,data,...){
  Roc.list(object=list(object),formula,data,...)
}

# }}}
# {{{ average Roc curves
avRoc <- function(list,grid,method="vertical"){
  if (missing(grid))
    grid <- switch(method,"vertical"=seq(0,1,.01),"horizontal"=seq(1,0,-.01))
  if (method=="vertical"){
    meanSens <- rowMeans(do.call("cbind",lapply(list,function(Roc){
      approx(x=Roc$Specificity,y=Roc$Sensitivity,xout=grid,ties=median,yleft=0,yright=1)$y
    })))
    meanRoc <- list(Sensitivity=c(meanSens,0),Specificity=c(grid,1))
  }
  else
    if (method=="horizontal"){
      meanSpec <- rowMeans(do.call("cbind",lapply(list,function(Roc){
        approx(x=Roc$Sensitivity,y=Roc$Specificity,xout=grid,ties=median,yleft=0,yright=1)$y
      })))
      meanRoc <- list(Sensitivity=c(grid,0),Specificity=c(meanSpec,1))
    }
  meanRoc
}
# }}}
# {{{ Roc.list

Roc.list <- function(object,
                     formula,
                     data,
                     splitMethod="noSplitMethod",
                     noinf.method=c("simulate"),
                     simulate="reeval",
                     B,
                     M,
                     breaks,
                     crRatio=1,
                     RocAverageMethod="vertical",
                     RocAverageGrid=switch(RocAverageMethod,
                         "vertical"=seq(0,1,.01),
                         "horizontal"=seq(1,0,-.01)),
                     model.args=NULL,
                     model.parms=NULL,
                     keepModels=FALSE,
                     keepSampleIndex=FALSE,
                     keepCrossValRes=FALSE,
                     keepNoInfSimu,
                     slaveseed,
                     na.accept=0,
                     verbose=TRUE,
                     ...){

    # }}}
    theCall=match.call()
    if (match("replan",names(theCall),nomatch=FALSE))
        stop("Argument name 'replan' has been replaced by 'splitMethod'.")
    
    # {{{ models
    NF <- length(object) 
  if (is.null(names(object)))names(object) <- sapply(object,function(o)class(o)[1])
  else{names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])}
  object.names = names(object)
  names(object) <- make.names(names(object),unique=TRUE)

# }}}
# {{{ formula
if (missing(formula)){
  formula <- eval(object[[1]]$call$formula)
  if (class(formula)!="formula")
    stop("Argument formula is missing.")
  else
    if (verbose)
      warning("Argument formula is missing. I use the formula from the call to the first model instead.")
}
  # }}}
# {{{ data
  if (missing(data)){
    data <- eval(object[[1]]$call$data)
    if (class(data)!="data.frame")
      stop("Argument data is missing.")
    else
      if (verbose)
        warning("Argument data is missing. I have (ab)used the data from the call\n of the first model instead.")
  
  }
  # }}}
# {{{ response

  m <- model.frame(formula,data,na.action=na.fail)
  Y <- model.response(m)
  
  if (is.factor(Y) && (length(levels(Y))==2) || length(unique(Y))==2) {
    Y <- factor(Y)
    Y <- as.numeric(Y==levels(Y)[2])
  }
  N <- length(Y)
  
    # }}}
    # {{{ break points for the ROC
    if (missing(breaks))
        if (splitMethod=="noSplitMethod")
            breaks <- NULL
        else
            breaks <- seq(0,1,.01)
    # }}}
  # {{{ SplitMethod
    SplitMethod <- MgSplitMethods(splitMethod=splitMethod,B=B,N=N,M=M,k=k)
    if (splitmethod$internal.name=="crossval")
        stop("K-fold cross-validation not available, but you can get similar results with bootstrap subsampling: set (1) splitMethod='bootcv' and (2) M = NROW(mydata)-round(0.1*NROW(mydata))")
    B <- SplitMethod$B
    CrossvalIndex <- SplitMethod$index
    if (!keepSampleIndex) SplitMethod$index <- NULL
    k <- SplitMethod$k
    do.crossval <- !(is.null(CrossvalIndex))
    if (missing(keepCrossValRes)) keepCrossValRes <- do.crossval
    if (missing(keepNoInfSimu)) keepNoInfSimu <- FALSE
    if (missing(slaveseed)||is.null(slaveseed))
        slaveseed <- sample(1:1000000,size=B,replace=FALSE)
    # }}}
    # {{{ checking the models for compatibility with cross-validation
    if (do.crossval){
        cm <- MgCheck(object=object,model.args=model.args,model.parms=model.parms,SplitMethod=SplitMethod,verbose=verbose)
        model.args <- cm$model.args
        model.parms <- cm$model.parms
    }
    # }}}
    # {{{ computation of ROC curves in a loop over the models 
    list.out <- lapply(1:NF,function(f){
        if (verbose && NF>1) message("\n",names(object)[f],"\n")
        fit <- object[[f]]
        extract <- model.parms[[f]]
    
  # }}}
# {{{ apparent ROC (use the same data for fitting and validation)

    pred <- do.call("predictStatusProb",c(list(object=fit,newdata=data),model.args[[f]]))
    AppRoc <- Roc.default(object=pred,y=Y,breaks=breaks,crRatio=crRatio)
    AppAuc <- Auc.default(object=AppRoc$Sensitivity,Spec=AppRoc$Specificity)
    AppBS <- Brier.default(object=pred,y=Y,crRatio=crRatio)
    
    # }}}
# {{{ No information error  
    if (SplitMethod$internal.name %in% c("boot632plus","noinf")){
        if (noinf.method=="simulate"){
            if (verbose)
                cat("\nSimulate no information performance\n")
            compute.NoInfRocList <- lapply(1:B,function(runb){
                if (verbose) MgTalk(runb,B)
                data.index <- data
                ## permute the response variable
                responseName <- all.vars(formula)[1]
                data.index[,responseName] <- sample(factor(Y),replace=FALSE)
                if (simulate=="reeval"){
                    fit.index <- MgRefit(object=fit,data=data.index,step=runb,silent=na.accept>0,verbose=verbose)
                }
                ## evaluate the model in data with reeallocated responses
                else
                    fit.index <- fit
                pred.index <- do.call("predictStatusProb",
                                      c(list(object=fit.index,newdata=data.index),
                                        model.args[[f]]))
                innerNoInfRoc <- Roc.default(object=pred.index,y=data.index[,responseName],breaks=breaks,crRatio=crRatio)
                innerNoInfBS <- Brier.default(object=pred.index,y=data.index[,responseName],crRatio=crRatio)
                list("innerNoInfRoc"=innerNoInfRoc,"innerNoInfBS"=innerNoInfBS)
            })
            if (verbose) cat("\n")
            NoInfRocList <- lapply(compute.NoInfRocList,function(x)x$innerNoInfRoc)
            NoInfRoc <- avRoc(list=NoInfRocList,grid=RocAverageGrid,method=RocAverageMethod)
            NoInfBS <- mean(sapply(compute.NoInfRocList,function(x){x$innerNoInfBS}))
            NoInfAuc <- mean(sapply(NoInfRocList,function(nil){Auc.default(object=nil$Sensitivity,nil$Specificity)}))
        }
        else{         
            NoInfRoc <- list(Sensitivity=c(breaks,0),Specificity=c(1-breaks,1))
            NoInfAuc <- 0.5
            NoInfBS <- .C("brier_noinf",bs=double(1),as.double(Y),as.double(pred),as.integer(N),NAOK=TRUE,PACKAGE="ModelGood")$bs
        }
    }
    if (SplitMethod$internal.name %in% c("boot632plus","bootcv","boot632")){
        # }}}
# {{{ Bootcv aka BootstrapCrossValidation
        if (verbose)
            cat("\nBootstrap cross-validation performance\n")

      compute.step <- function(runb,seed){
        if (verbose) MgTalk(runb,B)
        vindex.index <- match(1:N,CrossvalIndex[,runb],nomatch=0)==0
        val.index <- data[vindex.index,,drop=FALSE]
        train.index <- data[CrossvalIndex[,runb],,drop=FALSE]
        if (!is.null(seed)) {
          set.seed(seed)
        }
        fit.index <- MgRefit(object=fit,data=train.index,step=runb,silent=na.accept>0,verbose=verbose)
        if (!is.null(extract)) {
          fit.parms.index <- fit.index[extract]
          names(fit.parms.index) <- paste(extract,paste("sample",runb,sep="."),sep=":")
        }
        else fit.parms.index <- NULL
        if (is.null(fit.index)){
          failed <- "fit"
          innerBootcvRoc <- list(Sensitivity=NA,Specificity=NA)
          innerBCVBS <- NA
        }
        else{
          try2predict <- try(pred.index <- do.call("predictStatusProb",c(list(object=fit.index,newdata=val.index),model.args[[f]])),silent=na.accept>0)
          if (inherits(try2predict,"try-error")==TRUE){
            if (verbose) warning(paste("During bootstrapping: prediction for model ",class(fit.index)," failed in step ",runb),immediate.=TRUE)
            failed <- "prediction"
            innerBootcvRoc <- list(Sensitivity=NA,Specificity=NA)
            innerBCVBS <- NA
          }
          else{
            failed <- NA
            innerBootcvRoc <- Roc.default(y=Y[vindex.index],pred.index,breaks=breaks,crRatio=crRatio)
            innerBCVBS <- Brier.default(object=pred.index,y=Y[vindex.index],crRatio=crRatio)
          }
        }
        list("innerBootcvRoc"=innerBootcvRoc,
             "fit.parms"=fit.parms.index,
             "failed"=failed,
             "innerBCVBS"=innerBCVBS)
      }
      ## if (require(foreach)){
      b <- 1
      compute.BootcvRocList <- foreach (b = 1:B) %dopar% compute.step(runb=b,seed=slaveseed[[b]])
      ## }
      ## else{
      ## compute.BootcvRocList <- lapply(1:B,compute.step)
      ## }
      if (verbose) cat("\n")
      if (!is.null(extract)) fitParms <- sapply(compute.BootcvRocList,function(x)x$fit.parms)
      failed <- na.omit(sapply(compute.BootcvRocList,function(x)x$failed))
      BootcvRocList <- lapply(compute.BootcvRocList,function(x)x$innerBootcvRoc)
      BootcvRoc <- avRoc(BootcvRocList,grid=RocAverageGrid,method=RocAverageMethod)
      BCVSens <- BootcvRoc$Sensitivity
      BCVSpec <- BootcvRoc$Specificity
      BCVAuc <- mean(sapply(BootcvRocList,function(ool){Auc.default(object=ool$Sensitivity,ool$Specificity)}))
      BCVBS <- mean(sapply(compute.BootcvRocList,function(x){x$innerBCVBS}))
    }

# }}}
# {{{ Bootstrap .632

    if (SplitMethod$internal.name=="boot632"){
      B632Roc <- list(Sensitivity=.368 * AppRoc$Sensitivity + .632 * BCVSens,
                      Specificity=.368 * AppRoc$Specificity + .632 * BCVSpec)
      B632BS <- .368 * AppBS + .632 * BCVBS
      B632Auc <- .368 * AppAuc + .632 * BCVAuc
    }
  # }}}
# {{{ Bootstrap .632+
    if (SplitMethod$internal.name=="boot632plus"){
      ## first we have to prepare the averaging
      if (RocAverageMethod=="vertical"){
        AppSens <- c(approx(AppRoc$Sensitivity,AppRoc$Specificity,xout=RocAverageGrid,yleft=0,yright=1,ties=median)$y,0)
        AppRoc <- list(Sensitivity=AppSens,Specificity=c(RocAverageGrid,1))
        AppAuc <- Auc.default(object=AppRoc$Sensitivity,Spec=AppRoc$Specificity)
        
        R632Plus <- MgFormule632(App=AppSens,BCV=BCVSens,NoInf=NoInfRoc$Sensitivity,SmallerBetter=FALSE)
        B632plusRoc <- list(Sensitivity=R632Plus$B632Plus,Specificity=c(RocAverageGrid,1))
      }
      else{ ## RocAverageMethod horizontal
        AppSpec <- c(approx(AppRoc$Sensitivity,AppRoc$Specificity,xout=RocAverageGrid,yleft=1,yright=0,ties=median)$y,1)
        AppRoc <- list(Sensitivity=c(RocAverageGrid,0),Specificity=AppSpec)
        R632Plus <- MgFormule632(App=AppSpec,BCV=BCVSpec,NoInf=NoInfRoc$Specificity,SmallerBetter=FALSE)
        B632plusRoc <- list(Sensitivity=c(RocAverageGrid,0),
                            Specificity=R632Plus$B632Plus)
      }
      ## add the real .632+ AUC
      B632PlusAuc <- MgFormule632(App=AppAuc,BCV=BCVAuc,NoInf=NoInfAuc,SmallerBetter=FALSE)
      B632PlusBS <- MgFormule632(App=AppBS,BCV=BCVBS,NoInf=NoInfBS,SmallerBetter=TRUE)
    }
    # }}}
# {{{ preparing the output
    out <- switch(SplitMethod$internal.name,
                  "noSplitMethod"=list("Roc"=AppRoc),
                  ## "plain"=list("Roc"=BootRoc,"AppRoc"=AppRoc),
                  "boot632"=list("Roc"=B632Roc,"AppRoc"=AppRoc,"BootcvRoc"=BootcvRoc),
                  "boot632plus"=list("Roc"=B632plusRoc,
                    "AppRoc"=AppRoc,
                    "BootcvRoc"=BootcvRoc,
                    "NoInfRoc"=NoInfRoc,
                    "weight"=R632Plus$weight,
                    "overfit"=R632Plus$overfit),
                  "bootcv"=list("Roc"=BootcvRoc,"AppRoc"=AppRoc),
                  "noinf"=list("AppRoc"=AppRoc,"NoInfRoc"=NoInfRoc))
    
    out$Auc <- switch(SplitMethod$internal.name,
                      "noSplitMethod"=list("Auc"=AppAuc),
                      ## "plain"=list("Auc"=BootAuc,"AppAuc"=AppAuc),
                      "boot632"= list("Auc"=B632Auc,"AppAuc"=AppAuc,"AucBCV"=BCVAuc),
                      "boot632plus"=list("Auc"=B632PlusAuc$B632Plus,
                        "AppAuc"=AppAuc,
                        "BootcvAuc"=BCVAuc,
                        "NoInfAuc"=NoInfAuc,
                        "weight"=B632PlusAuc$weight,
                        "overfit"=B632PlusAuc$overfit),
                      "bootcv"=list("Auc"=BCVAuc,"AppAuc"=AppAuc),
                      "noinf"=list("NoInfAuc"=NoInfAuc,"AppAuc"=AppAuc))
    out$Brier <- switch(SplitMethod$internal.name,
                        "noSplitMethod"=list("BS"=AppBS),
                        ## "plain"=list("BS"=BootBS,"AppBS"=AppBS),
                        "boot632"= list("BS"=B632BS,"AppBS"=AppBS,"BSBCV"=BCVBS),
                        "boot632plus"=list("BS"=B632PlusBS$B632Plus,
                          "AppBS"=AppBS,
                          "BootcvBS"=BCVBS,
                          "NoInfBS"=NoInfBS,
                          "weight"=B632PlusBS$weight,
                          "overfit"=B632PlusBS$overfit),
                        "bootcv"=list("BS"=BCVBS,"AppBS"=AppBS),
                        "noinf"=list("AppBS"=AppBS,"NoInfBS"=NoInfBS))

    ##     if (keepCrossValRes==TRUE && SplitMethod$internal.name!="noSplitMethod"){
    
    if (keepCrossValRes==TRUE && class(try(is.null(BootcvRocList),silent=TRUE))!="try-error"){
      if (SplitMethod$internal.name!="noinf")
        out <- c(out,list("BootcvRocList"=BootcvRocList))
    }
    if (keepNoInfSimu==TRUE && SplitMethod$internal.name!="noSplitMethod"){
      if (SplitMethod$internal.name!="noinf")
        out <- c(out,list("NoInfRocList"=NoInfRocList))
    }
    if (!is.null(extract)) out <- c(out,list("fitParms"=fitParms))
    if (na.accept>0) out <- c(out,list("failed"=failed))
    out
  })
  ## it might be that the first model has no extracted parameters
  ## but one of the other has
  if (length(model.parms)>0)
    names.lout <- c(names(list.out[[1]]),"fitParms")
  else
    names.lout <- names(list.out[[1]])
  out <- lapply(names.lout,function(w){
    e <- lapply(list.out,function(x){x[[w]]})
    names(e) <- names(object)
    e
  })
  names(out) <- names.lout
  if(keepModels==TRUE)
    outmodels <- object
  else if (keepModels=="Call"){
    outmodels <- lapply(object,function(o)o$call)
    names(outmodels) <- names(object)
  }
  else{
    outmodels <- names(object)
    names(outmodels) <- names(object)
  }
  out <- c(out,
           list(call=match.call(),
                Response=Y,
                models=outmodels,
                model.names=object.names,
                method=SplitMethod,
                breaks=breaks,crRatio=crRatio))
  if (verbose) cat("\n")
  # }}}
class(out) <- "Roc"
out
}

