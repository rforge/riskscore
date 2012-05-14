score <- function(object,...){
  UseMethod("score",object=object)
}
score.list <- function(object,
                       formula,
                       data,
                       refLevel,
                       measures=c("Auc","Brier","cindex","Roc"),
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
                       ...){
  
  theCall <- match.call(expand=TRUE)
  # ------------------parse arguments and prepare data------------------
  # {{{ Response
  formula <- getFormula(formula,object,verbose=verbose)
  data <- getData(data=data,formula=formula,object=object,verbose=verbose)
  response <- getResponse(formula=update(formula,~1),data=data,refLevel=refLevel,verbose=verbose)
  N <- NROW(response)
  responseType <- attr(response,"model")
  predictHandlerFun <- switch(responseType,"survival"="predictSurvProb","competing.risks"="predictEventProb","binary"="predictStatusProb",stop("Dont know how to predict response of type ",responseType))
  censType <- attr(response,"cens.type")
  if (is.null(censType)) censType <- "uncensoredData"
  # }}}
  # {{{ SplitMethod
  splitMethod <- getSplitMethod(splitMethod=splitMethod,B=B,N=N,M=M)
  splitIndex <- splitMethod$index
  do.resample <- !(is.null(splitIndex))
  # }}}
  # {{{ Prediction models
  # checking the models for predictHandlerFunction
  allmethods <- methods(predictHandlerFun)
  wantedMethods <- paste(predictHandlerFun,sapply(object,class),sep=".")
  if (any(notfound <- match(wantedMethods,allmethods,nomatch=FALSE)==0))
    stop(paste("Could not find all predict methods for",responseType,"response:",paste(wantedMethods[notfound],sep=", ")))
  
  # checking the models for compatibility with resampling
  if (do.resample){
    object <- lapply(1:length(object),function(f){
      fit <- object[[f]]
      if(is.null(fit$call))
        stop(paste("model",names(object)[f],"does not have a call argument."))
      else fit$call$data <- NULL
    })
  }
  # add null model and find names for the object
  if (nullModel==TRUE)
    object <- c(getNullModel(formula=formula,data=data,responseType=responseType),object)
  if (is.null(names(object)))
    {names(object) <- sapply(object,function(o)class(o)[1])}
  else
    {names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
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
  # -----------------deal with censored data---------------------------
  # {{{ Censoring
  if (censType=="rightCensored"){
    if (censMethod=="ipcw"){
      weights <- getCensoringWeights(formula=formula,
                                     data=data,
                                     times=times,
                                     censModel=censModel,
                                     responseType=responseType)
    }
    else{
      censMethod <- "jackknife.pseudo.values"
      pseudoResponse <- getPseudoValues(formula=formula,
                                        data=data,
                                        responseType=responseType,
                                        times=times,
                                        cause=cause)
    }
  }
  else{
    censMethod <- "none"
  }
  # }}}
  # ---------------------to come to the conclusion---------------------
  # {{{ Handler functions
  if (verbose==1){
    cat(paste("\nresponseType:",responseType))
    cat(paste("\npredictHandlerFun:",predictHandlerFun))
    cat(paste("\ncensType:",censType))
    cat(paste("\ncensoringHandling:",censMethod))
    cat(paste("\nsplitMethod:",splitMethod$name))
    cat("\n")
  }
  # }}}
  # ------------------------------metrics------------------------------
  # {{{
  useMeasures <- na.omit(match(tolower(measures),tolower(c("AUC","Brier","Roc","cindex"))))
  # }}}
  # -----------------apparent performance-----------------------
  # {{{  
  browser()
  # }}}
  #------------------output------------------------------------
  
}


