##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param object 
##' @param ... 
##' @return 
##' @author Thomas Alexander Gerds
score <- function(object,...){
  UseMethod("score",object=object)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param object 
##' @param formula 
##' @param data 
##' @param event 
##' @param measures 
##' @param times 
##' @param maxtime 
##' @param landmarks 
##' @param useObservedTimes 
##' @param nullModel 
##' @param censMethod 
##' @param censModel 
##' @param splitMethod 
##' @param B 
##' @param M 
##' @param verbose 
##' @param ... 
##' @return 
##' @author Thomas Alexander Gerds
score.list <- function(object,
                       formula,
                       data,
                       event,
                       measures=c("Auc","Brier","cindex","Roc","calPlot"),
                       times,
                       maxtime,
                       landmarks,
                       useObservedTimes=TRUE,
                       nullModel=TRUE,
                       censMethod="ipcw",
                       censModel="cox",
                       splitMethod,
                       B,
                       M,
                       verbose=1,
                       trainseeds,
                       ...){
  
  theCall <- match.call(expand=TRUE)
  # -----------------parse arguments and prepare data---------
  # {{{ Response
  formula <- getFormula(formula,object,verbose=verbose)
  data <- getData(data=data,formula=formula,object=object,verbose=verbose)
  responseFormula <- update(formula,~1)
  if (missing(event)) event <- NULL
  response <- getResponse(formula=responseFormula,data=data,event=event,verbose=verbose)
  N <- NROW(response)
  responseType <- attr(response,"model")
  predictHandlerFun <- switch(responseType,"survival"="predictSurvProb","competing.risks"="predictEventProb","binary"="predictStatusProb",stop("Dont know how to predict response of type ",responseType))
  censType <- attr(response,"cens.type")
  if (is.null(censType)) censType <- "uncensoredData"
  # }}}
  # {{{ SplitMethod
  splitMethod <- getSplitMethod(splitMethod=splitMethod,B=B,N=N,M=M)
  B <- splitMethod$B
  splitIndex <- splitMethod$index
  do.resample <- !(is.null(splitIndex))
  # }}}
  # {{{ Checking the models for predictHandlerFunction
  allmethods <- methods(predictHandlerFun)
  wantedMethods <- lapply(object,function(o){
    candidateMethods <- paste(predictHandlerFun,class(o),sep=".")
    if (all(match(candidateMethods,allmethods,nomatch=0)==0))
    stop(paste("Could not find ",predictHandlerFun," method for ",paste(class(o),collapse=" ,"),sep=""))
  })
  # checking the models for compatibility with resampling
  if (do.resample){
    object <- lapply(1:length(object),function(f){
      fit <- object[[f]]
      if(is.null(fit$call))
        stop(paste("model",names(object)[f],"does not have a call argument."))
      else fit$call$data <- NULL
      fit
    })
  }
  # add null model and find names for the object
  if (nullModel==TRUE)
    object <- c(getNullModel(formula=formula,data=data,responseType=responseType),object)
  if (is.null(names(object))){
      names(object) <- sapply(object,function(o)class(o)[1])}
  else {
    names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
   }
  names(object) <- make.names(names(object),unique=TRUE)
  NF <- length(object)

  # }}}
  # {{{ Evaluation landmarks and horizons (times)

  if (responseType %in% c("survival","competing.risks")){
    if (missing(maxtime) || is.null(maxtime)){
      eventTimes <- unique(response[,"time"])
      maxtime <- eventTimes[length(eventTimes)]
    }
    if (missing(landmarks)){
      start <- 0
      if (missing(times)){
        if (useObservedTimes==TRUE)
          times <- unique(c(start,eventTimes))
        else
          times <- seq(start,maxtime,(maxtime - start)/100)
      }
      else{
        if (useObservedTimes==TRUE) 
          times <- sort(c(start,unique(times),eventTimes))
        else
          times <- sort(unique(c(start,times)))
      }
      times <- times[times<=maxtime]
      NT <-  length(times)
    }
    else{
      stop("Landmark updating not yet implemented.")
    }
  }
  else{
    if (!missing(times)) warning("Function 'score': Response type is not time-to-event: argument 'times' will be ignored.",call.=FALSE)
  }

  # }}}
  # -----------------deal with censored data------------------
  # {{{ 

  if (censType=="rightCensored"){
    if (censMethod=="ipcw"){
      weights <- getCensoringWeights(formula=formula,data=data,times=times,censModel=censModel,responseType=responseType)
    }
    else{
      censMethod <- "jackknife.pseudo.values"
      pseudoResponse <- getPseudoValues(formula=formula,data=data,responseType=responseType,times=times,event=event)
    }
  }
  else{
    censMethod <- "none"
  }

  # }}}
  # -----------------specific task handlers-------------------
  # {{{ 
  if (verbose==1){
    cat(paste("\nresponseType:",responseType))
    cat(paste("\npredictHandlerFun:",predictHandlerFun))
    cat(paste("\ncensType:",censType))
    cat(paste("\ncensoringHandling:",censMethod))
    cat(paste("\nsplitMethod:",splitMethod$name))
    cat("\n")
  }
  # prediction list
  plist <- switch(predictHandlerFun,"predictStatusProb"={list(newdata=data)},
                  "predictSurvProb"={list(newdata=data,times=times)},
                  "predictEventProb"={list(newdata=data,times=times,event=event)},
                  stop("Unknown predictHandler."))
  # }}}
  # -----------------output structure: performance----------------------
  # {{{

  ## Metrics <- na.omit(match(tolower(measures),tolower(c("AUC","Brier","cindex","Roc","calPlot"))))
  Metrics <- lapply(measures,grep,c("AUC","Brier","cindex","Roc","calPlot"),ignore.case=TRUE,value=TRUE)
  performance <- lapply(Metrics,function(m){
    vector(NF,mode="list")
  })
  names(performance) <- Metrics
  # metrics list
  mlist <- list(response=response,"prediction"=NULL)
  if (predictHandlerFun %in% c("predictSurvProb","predictEventProb")){
    mlist <- c(mlist,list(times=times))
    if (censMethod=="ipcw"){
      mlist <- c(mlist,list(weights=weights))
    }
  }

  # }}}
  # -----------------apparent performance---------------------
  # {{{
  calcPerf <- function(b=0,traindata,trainseed,testdata){
    plist[["newdata"]] <- testdata
    model.perf <- lapply(1:NF,function(f){
      if (!is.null(traindata)){
        set.seed(trainseed)
        trained.model <- trainModel(model=object[[f]],data=traindata)
      }
      else{
        trained.model <- object[[f]]
      }
      pred <- do.call(predictHandlerFun,c(list(object=trained.model),plist))
      mlist[["prediction"]] <- pred
      mlist[["response"]] <- getResponse(formula=responseFormula,data=testdata,event=event,verbose=verbose)
      metric <- lapply(c("Brier","AUC"),function(m){
        do.call(paste(m,responseType,sep="."),mlist)
      })
      names(metric) <- "Brier"
      metric
    })
    names(model.perf) <- names(object)
    model.perf
  }
  trainModel <- function(model,data){
    model$call$data <- data
    try(eval(model$call),silent=TRUE)
  }
  noSplit <- calcPerf(b=0,traindata=NULL,testdata=data)
  crossval <- NULL
  if (splitMethod$name=="BootCv"){
    message("Running cross-validation")
    if (missing(trainseeds)||is.null(trainseeds))
      trainseeds <- sample(1:1000000,size=B,replace=FALSE)
    if (require(foreach)){
      crossval <- foreach (b=1:B) %dopar% {
        calcPerf(b,
                 traindata=data[splitMethod$index[,b],,drop=FALSE],
                 trainseed=trainseeds[b],
                 testdata=data[match(1:N,unique(splitMethod$index[,b]),nomatch=0)==0,,drop=FALSE])
      }
    }
  }
  # }}}
  #------------------output-----------------------------------
  output <- list(noSplitPerf=noSplit,crossValPerf=crossval)
}


Brier.binary <- function(response,prediction){
  mean((response-prediction)^2)
}

AUC.binary <- function(response,prediction,breaks,costs=1,yes=FALSE){
  N <- length(prediction)
  if(length(response)!=N) stop("Arguments must have the same length")
  if(length(unique(response))!=2) stop("response must be binary")
  Event <- as.integer(as.character(factor(response,labels=c("0","1"))))
  count.EventPos <- sum(Event==1)
  count.EventNeg <- sum(Event==0)
  if (missing(breaks))
    breaks <- sort(unique(prediction))
  else
    breaks <- sort(unique(breaks))
  if (length(breaks)>1 & !is.factor(prediction)){
    eval.times <- breaks- min(diff(breaks))/2
    count.TestPos <- N-sindex(jump.times=prediction,eval.times=eval.times)
    count.TestNeg <- N-count.TestPos
    a <- count.EventPos-sindex(jump.times=prediction[Event==1],eval.times=eval.times)
    b <- count.EventNeg-sindex(jump.times=prediction[Event==0],eval.times=eval.times)
  }
  else{
    tabx <- table(prediction)
    count.TestPos <- tabx
    count.TestNeg <- tabx
    tabxpos <- table(prediction[Event==1])
    tabxneg <- table(prediction[Event==0])
    a <- count.EventPos-tabxpos
    b <- count.EventNeg-tabxneg
  }
  c <- count.EventPos-a
  d <- count.EventNeg-b
  Sens <- c(costs*a/count.EventPos,0)
  Spec <- c(d/count.EventNeg,1)
  N <- length(Sens)
  if (length(unique(Spec)) < 2) return(NA)
  if (Spec[1]<Spec[N]) {Spec <- rev(Spec); Sens <- rev(Sens)}
  auc <- 0.5 * sum(diff(c(0,1-Spec,1)) * (c(Sens,1) + c(0,Sens)))
  auc
}
