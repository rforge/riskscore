CindexBootstrapCrossValidation <- function(object,
                                           data,
                                           Y,
                                           status,
                                           event,
                                           eval.times,
                                           pred.times,
                                           cause,
                                           weights,
                                           ipcw.refit=FALSE,
                                           ipcw.call,
                                           tiedPredictionsIn,
                                           tiedOutcomeIn,
                                           tiedMatchIn,
                                           splitMethod,
                                           multiSplitTest,
                                           keepResiduals,
                                           testTimes,
                                           confInt,
                                           confLevel,
                                           getFromModel,
                                           giveToModel,
                                           predictHandlerFun,
                                           keepMatrix,
                                           verbose,
                                           savePath,slaveseed){
  
  # {{{ initializing
  B <- splitMethod$B
  N <- splitMethod$N
  M <- splitMethod$M
  NT <- length(eval.times)
  NF <- length(object) 
  ResampleIndex <- splitMethod$index
  # }}}
  step <- function(b,seed){
    if (verbose==TRUE) internalTalk(b,B)
    # {{{ training and validation data
    vindex.b <- match(1:N,unique(ResampleIndex[,b]),nomatch=0)==0
    Y.b <- Y[vindex.b]
    tindex.b <- match(Y.b,unique(Y.b))
    val.b <- data[vindex.b,,drop=FALSE]
    ## browser()
    train.b <- data[ResampleIndex[,b],,drop=FALSE]
    ## if (b==1) print(train.b$days)
    ## if (b==1) print(val.b$days)
    NV=sum(vindex.b)                    # NROW(val.b)
    # }}}
    # {{{ IPCW
    ## if (ipcw.refit==TRUE){
    ## ipcw.call.b.i <- ipcw.call$weight.i
    ## ipcw.call.b.j <- ipcw.call$weight.j
    ## ipcw.call.b.i$data <- val.b
    ## ipcw.call.b.j$data <- val.b
    ## ipcw.call.b.i$subjectTimes <- Y.b
    ## ipcw.call.b.j$subjectTimes <- Y.b
    ## ipcw.b.i <- do.call("ipcw",ipcw.call.b.i)$IPCW.subjectTimes
    ## ipcw.b.j <- do.call("ipcw",ipcw.call.b.j)$IPCW.times
    ## }
    ## else{
    ipcw.b.i <- weights$weight.i[vindex.b]
    if (is.null(dim(weights$weight.j))){
      ipcw.b.j <- weights$weight.j
    }
    else{
      ipcw.b.j <- weights$weight.j[vindex.b,]
    }
    ## }
    # }}}
    # {{{ Building the models in training data
    if (!is.null(seed)) {
      set.seed(seed)
      ## message("seed:",seed)
    }
    trainModels <- lapply(1:NF,function(f){
      fit.b <- internalReevalFit(object=object[[f]],
                                 data=train.b,
                                 step=b,
                                 silent=FALSE,
                                 verbose=verbose)
      ## fit.b$call <- object[[f]]$call
      fit.b
    })
    # }}}
    # {{{ Saving the models?
    if (!is.null(savePath)){
      nix <- lapply(1:NF,function(f){
        fit.b <- trainModels[[f]]
        ## print(object.size(fit.b))
        fit.b$formula <- NULL
        ## print(environment(fit.b$formula))
        save(fit.b,file=paste(paste(savePath,"/",names(object)[f],"-bootstrap-",b,sep=""),".rda",sep=""))
      })
    }
    # }}}
    # {{{ Extracting parameters?
    if (!is.null(getFromModel)){
      ModelParameters <- lapply(1:NF,function(f){
        getParms <- getFromModel[[f]]
        if (is.null(getParms)) trainModels[[f]][getParms] else NULL
      })
    }
    # }}}
    # {{{ Check fits
    fitFailed <- lapply(trainModels,function(fit.b) (is.null(fit.b)))
    # }}}
    # {{{ Predicting the validation data
    predVal <- lapply(1:NF,function(f){
      fit.b <- trainModels[[f]]
      extraArgs <- giveToModel[[f]]
      if (predictHandlerFun %in% c("predictEventProb","predictLifeYearsLost")){
        try2predict <- try(pred.b <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=val.b,times=pred.times,train.data=train.b,cause=cause),extraArgs)))
      }
      else{
        try2predict <- try(pred.b <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=val.b,times=pred.times,train.data=train.b),extraArgs)))
      }
      ## browser()
      ## print(pred.b[1:5])
      if (inherits(try2predict,"try-error")==TRUE){
        if (verbose==TRUE) warning(paste("During bootstrapping: prediction for model ",class(fit.b)," failed in step ",b),immediate.=TRUE)
        NULL}
      else{
        pred.b
      }
    })
    # }}}
    # {{{ Compute cindex  for step b
    if (multiSplitTest==TRUE){
      stop("not yet defined: residual test for cindex")
      Residuals <- lapply(predVal,function(pred.b){
        if (is.null(pred.b))
          NA
        else{
          if (predictHandlerFun %in% c("predictEventProb","predictLifeYearsLost")){
            1
            ## matrix(.C("pecResidualsCR",pec=double(NT),resid=double(NT*NV),as.double(Y[vindex.b]),as.double(status[vindex.b]),as.double(event[vindex.b]),as.double(times),as.double(pred.b),as.double(ipcwTimes.b),as.double(IPCW.subjectTimes.b),as.integer(NV),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred.b))),NAOK=TRUE,PACKAGE="pec")$resid,ncol=NT,byrow=FALSE)
          }
          else{
            1
            ## matrix(.C("pecResiduals",pec=double(NT),resid=double(NT*NV),as.double(Y[vindex.b]),as.double(status[vindex.b]),as.double(times),as.double(pred.b),as.double(ipcwTimes.b),as.double(IPCW.subjectTimes.b),as.integer(NV),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred.b))),NAOK=TRUE,PACKAGE="pec")$resid,ncol=NT,byrow=FALSE)
          }
        }
      })
      names(Residuals) <- names(object)
      ## PredCindexStepB=lapply(Residuals,function(x){colMeans(x)})
      PredCindexStepB=1
    }
    else{
      PredCindexStepB <- lapply(predVal,function(pred.b){
        if (is.null(pred.b))
          NA
        else{
          if (predictHandlerFun %in% c("predictEventProb","predictLifeYearsLost")){
            Step.b.CindexResult <- .C("ccr",cindex=double(NT),concA=double(NT),pairsA=double(NT),concB=double(NT),pairsB=double(NT),as.integer(tindex.b),as.double(Y.b),as.integer(status[vindex.b]),as.integer(event[vindex.b]),as.double(eval.times),as.double(ipcw.b.i),as.double(ipcw.b.j),as.double(pred.b),as.integer(sum(vindex.b)),as.integer(NT),as.integer(tiedPredictionsIn),as.integer(tiedOutcomeIn),as.integer(tiedMatchIn),as.integer(!is.null(dim(ipcw.b.j))),NAOK=TRUE,package="pec")
            Step.b.Cindex <- Step.b.CindexResult$cindex
            Step.b.PairsA <- Step.b.CindexResult$pairsA
            Step.b.ConcordantA <- Step.b.CindexResult$concA
            Step.b.PairsB <- Step.b.CindexResult$pairsB
            Step.b.ConcordantB <- Step.b.CindexResult$concB
            list(Cindex.b=Step.b.Cindex,Pairs.b=list(A=Step.b.PairsA,B=Step.b.PairsB),Concordant.b=list(A=Step.b.ConcordantA,B=Step.b.ConcordantB))
          }
          else{
            cindexOut <- .C("cindex",


#' Concordance index for right censored survival time data
#' 
#' In survival analysis, a pair of patients is called concordant if the risk of
#' the event predicted by a model is lower for the patient who experiences the
#' event at a later timepoint. The concordance probability (C-index) is the
#' frequency of concordant pairs among all pairs of subjects. It can be used to
#' measure and compare the discriminative power of a risk prediction models.
#' The function provides an inverse of the probability of censoring weigthed
#' estimate of the concordance probability to adjust for right censoring.
#' Cross-validation based on bootstrap resampling or bootstrap subsampling can
#' be applied to assess and compare the discriminative power of various
#' regression modelling strategies on the same set of data.
#' 
#' Pairs with identical observed times, where one is uncensored and one is
#' censored, are always considered usuable (independent of the value of
#' \code{tiedOutcomeIn}), as it can be assumed that the event occurs at a later
#' timepoint for the censored observation.
#' 
#' For uncensored response the result equals the one obtained with the
#' functions \code{rcorr.cens} and \code{rcorrcens} from the \code{Hmisc}
#' package (see examples).
#' 
#' @aliases cindex cindex.list
#' @param object A named list of prediction models, where allowed entries are
#' (1) R-objects for which a \link{predictSurvProb} method exists (see
#' details), (2) a \code{call} that evaluates to such an R-object (see
#' examples), (3) a matrix with predicted probabilities having as many rows as
#' \code{data} and as many columns as \code{times}. For cross-validation all
#' objects in this list must include their \code{call}.
#' @param formula A survival formula. The left hand side is used to finde the
#' status response variable in \code{data}. For right censored data, the right
#' hand side of the formula is used to specify conditional censoring models.
#' For example, set \code{Surv(time,status)~x1+x2} and \code{cens.model="cox"}.
#' Then the weights are based on a Cox regression model for the censoring times
#' with predictors x1 and x2.  Note that the usual coding is assumed:
#' \code{status=0} for censored times and that each variable name that appears
#' in \code{formula} must be the column name in \code{data}. If there are no
#' covariates, i.e. \code{formula=Surv(time,status)~1} the \code{cens.model} is
#' coerced to \code{"marginal"} and the Kaplan-Meier estimator for the
#' censoring times is used to calculate the weights.  If \code{formula} is
#' \code{missing}, try to extract a formula from the first element in object.
#' @param data A data frame in which to validate the prediction models and to
#' fit the censoring model.  If \code{data} is missing, try to extract a data
#' set from the first element in object.
#' @param eval.times A vector of timepoints for evaluating the discriminative
#' ability. At each timepoint the c-index is computed using only those pairs
#' where one of the event times is known to be earlier than this timepoint. If
#' \code{eval.times} is \code{missing} or \code{Inf} then the largest
#' uncensored event time is used.
#' @param pred.times A vector of timepoints for evaluating the prediction
#' models. This should either be exactly one timepoint used for all
#' \code{eval.times}, or be as long as \code{eval.times}, in which case the
#' predicted order of risk for the jth entry of \code{eval.times} is based on
#' the jth entry of \code{pred.times} corresponding
#' @param cause For competing risks, the event of interest. Defaults to the
#' first state of the response, which is obtained by evaluating the left hand
#' side of \code{formula} in \code{data}.
#' @param lyl If \code{TRUE} rank subjects accoring to predicted
#' life-years-lost (See Andersen due to this cause instead of predicted risk.
#' @param cens.model Method for estimating inverse probability of censoring
#' weigths:
#' 
#' \code{cox}: A semi-parametric Cox proportional hazard model is fitted to the
#' censoring times
#' 
#' \code{marginal}: The Kaplan-Meier estimator for the censoring times
#' 
#' \code{nonpar}: Nonparametric extension of the Kaplan-Meier for the censoring
#' times using symmetric nearest neighborhoods -- available for arbitrary many
#' strata variables on the right hand side of argument \code{formula} but at
#' most one continuous variable. See the documentation of the functions
#' \code{prodlim} and \code{neighborhood} from the prodlim package.
#' 
#' \code{aalen}: The nonparametric Aalen additive model fitted to the censoring
#' times. Requires the timereg package maintained by Thomas Scheike.
#' @param ipcw.refit If \code{TRUE} the inverse probability of censoring
#' weigths are estimated separately in each training set during
#' cross-validation.
#' @param ipcw.limit Value between 0 and 1 (but no equal to 0!) used to cut for
#' small weights in order to stabilize the estimate at late times were few
#' individuals are observed.
#' @param tiedPredictionsIn If \code{FALSE} pairs with identical predictions
#' are excluded, unless also the event times are identical and uncensored and
#' \code{tiedMatchIn} is set to \code{TRUE}.
#' @param tiedOutcomeIn If \code{TRUE} pairs with identical and uncensored
#' event times are excluded, unless also the predictions are identical and
#' \code{tiedMatchIn} is set to \code{TRUE}.
#' @param tiedMatchIn If \code{TRUE} then pairs with identical predictions and
#' identical and uncensored event times are counted as concordant pairs.
#' @param splitMethod Defines the internal validation design:
#' 
#' \code{none/noPlan}: Assess the models in the give \code{data}, usually
#' either in the same data where they are fitted, or in independent test data.
#' 
#' \code{BootCv}: Bootstrap cross validation. The prediction models are trained
#' on \code{B} bootstrap samples, that are either drawn with replacement of the
#' same size as the original data or without replacement from \code{data} of
#' the size \code{M}.  The models are assessed in the observations that are NOT
#' in the bootstrap sample.
#' 
#' \code{Boot632}: Linear combination of AppCindex and OutOfBagCindex using the
#' constant weight .632.
#' 
#' @param B Number of bootstrap samples. The default depends on argument
#' \code{splitMethod}.  When \code{splitMethod} in c("BootCv","Boot632") the
#' default is 100.  For \code{splitMethod="none"} \code{B} is the number of
#' bootstrap simulations e.g. to obtain bootstrap confidence limits -- default
#' is 0.
#' @param M The size of the bootstrap samples for resampling without
#' replacement. Ignored for resampling with replacement.
#' @param model.args List of extra arguments that can be passed to the
#' \code{predictSurvProb} methods. The list must have an entry for each entry
#' in \code{object}.
#' @param model.parms Experimental.  List of with exactly one entry for each
#' entry in \code{object}.  Each entry names parts of the value of the fitted
#' models that should be extracted and added to the value.
#' @param keep.index Logical. If \code{FALSE} remove the bootstrap or
#' cross-validation index from the output list which otherwise is included in
#' the method part of the output list.
#' @param keep.matrix Logical. If \code{TRUE} add all \code{B} prediction error
#' curves from bootstrapping or cross-validation to the output.
#' @param keep.models Logical. If \code{TRUE} keep the models in object. If
#' \code{"Call"} keep only the \code{call} of these models.
#' @param keep.residuals Experimental.
#' @param keep.pvalues Experimental.
#' @param keep.weights Experimental.
#' @param multiSplitTest Experimental.
#' @param testTimes A vector of time points for testing differences between
#' models in the time-point specific Brier scores.
#' @param confInt Experimental.
#' @param confLevel Experimental.
#' @param verbose if \code{TRUE} report details of the progress, e.g. count the
#' steps in cross-validation.
#' @param savePath Place in your filesystem (directory) where training models
#' fitted during cross-validation are saved. If \code{missing} training models
#' are not saved.
#' @param slaveseed Vector of seeds, as long as \code{B}, to be given to the
#' slaves in parallel computing.
#' @param na.action Passed immediately to model.frame. Defaults to na.fail. If
#' set otherwise most prediction models will not work.
#' @param ... Not used.
#' @return Estimates of the C-index.
#' @author Thomas A Gerds \email{tag@@biostat.ku.dk}
#' @references
#' 
#' TA Gerds, MW Kattan, M Schumacher, and C Yu. Estimating a time-dependent
#' concordance index for survival prediction models with covariate dependent
#' censoring. Statistics in Medicine, Ahead of print:to appear, 2013. DOI =
#' 10.1002/sim.5681
#' 
#' Wolbers, M and Koller, MT and Witteman, JCM and Gerds, TA (2013) Concordance
#' for prognostic models with competing risks Research report 13/3. Department
#' of Biostatistics, University of Copenhagen
#' 
#' Andersen, PK (2012) A note on the decomposition of number of life years lost
#' according to causes of death Research report 12/2. Department of
#' Biostatistics, University of Copenhagen
#' @keywords survival
#' @examples
#' 
#' \donttest{
#'  # simulate data based on Weibull regression  
#' library(prodlim)
#'  set.seed(13)
#'  dat <- SimSurv(300)
#'  # fit three different Cox models and a random survival forest
#'  # note: low number of trees for the purpose of illustration 
#'  library(survival)
#'  library(randomForestSRC)
#'  cox12 <- coxph(Surv(time,status)~X1+X2,data=dat)
#'  cox1 <- coxph(Surv(time,status)~X1,data=dat)
#'  cox2 <- coxph(Surv(time,status)~X2,data=dat)
#'  rsf1 <- rfsrc(Surv(time,status)~X1+X2,data=dat,ntree=15,forest=TRUE)
#'  #
#'  # compute the apparent estimate of the C-index at different time points
#'  #
#' ApparrentCindex  <- cindex(list("Cox X1"=cox1,
#' 		       "Cox X2"=cox2,
#' 		       "Cox X1+X2"=cox12,
#'                        "RSF"=rsf1),
#' 		  formula=Surv(time,status)~X1+X2,
#' 		  data=dat,
#' 		  eval.times=seq(5,500,50))
#'   print(ApparrentCindex)
#'   plot(ApparrentCindex)
#'  #
#'  # compute the bootstrap-crossvalidation estimate of
#'  # the C-index at different time points
#'  #
#' set.seed(142)
#' bcvCindex  <- cindex(list("Cox X1"=cox1,
#' 		       "Cox X2"=cox2,
#' 		       "Cox X1+X2"=cox12,
#'                        "RSF"=rsf1),
#' 		  formula=Surv(time,status)~X1+X2,
#' 		  data=dat,
#'                   splitMethod="bootcv",
#'                   B=10,
#' 		  eval.times=seq(5,500,50))
#'   print(bcvCindex)
#'   plot(bcvCindex)
#' }
#' \donttest{
#'  # for uncensored data the results are the same
#'  # as those obtained with the function rcorr.cens from Hmisc
#' library(Hmisc)
#' set.seed(16)
#' dat <- SimSurv(30,cens=FALSE)
#' fit12 <- coxph(Surv(time,status)~X1+X2,data=dat)
#' fit1 <- coxph(Surv(time,status)~X1,data=dat)
#' fit2 <- coxph(Surv(time,status)~X2,data=dat)
#' Cpec <- cindex(list("Cox X1+X2"=fit12,"Cox X1"=fit1,"Cox X2"=fit2),
#' 	       formula=Surv(time,status)~1,
#' 	       data=dat,
#' 	       eval.times=Inf)
#' p1 <- predictSurvProb(fit1,newdata=dat,times=100)
#' p2 <- predictSurvProb(fit2,newdata=dat,times=100)
#' p12 <- predictSurvProb(fit12,newdata=dat,times=100)
#' harrelC1 <- rcorr.cens(p1,with(dat,Surv(time,status)))
#' harrelC2 <- rcorr.cens(p2,with(dat,Surv(time,status)))
#' harrelC12 <- rcorr.cens(p12,with(dat,Surv(time,status)))
#' harrelC1[["C Index"]]==Cpec$AppCindex[["Cox.X1"]]
#' harrelC2[["C Index"]]==Cpec$AppCindex[["Cox.X2"]]
#' harrelC12[["C Index"]]==Cpec$AppCindex[["Cox.X1.X2"]]
#'  #
#'  # competing risks 
#'  #
#' # require(lava)
#' % m <- lvm()
#' % X <- paste("X",seq(5),sep="")
#' 
#' }
#' 
#' @export cindex
                            cindex=double(NT),
                            conc=double(NT),
                            pairs=double(NT),
                            as.integer(tindex.b),
                            as.double(Y.b),
                            as.integer(status[vindex.b]),
                            as.double(eval.times),
                            as.double(ipcw.b.i),
                            as.double(ipcw.b.j),
                            as.double(pred.b),
                            as.integer(sum(vindex.b)),
                            as.integer(NT),
                            as.integer(tiedPredictionsIn),
                            as.integer(tiedOutcomeIn),
                            as.integer(tiedMatchIn),
                            as.integer(!is.null(dim(ipcw.b.j))),
                            NAOK=TRUE,
                            package="pec")            
            Cindex.b <- cindexOut$cindex
            Pairs.b <- cindexOut$pairs 
            Concordant.b <- cindexOut$conc
            list(Cindex.b=Cindex.b,Pairs.b=Pairs.b,Concordant.b=Concordant.b)
          }
        }
      })
    }
    # }}}
    # {{{ van de Wiel's test
    ## if (multiSplitTest==TRUE){
    ## testedResid <- testResiduals(Residuals,times=times,testTimes=testTimes,rangeInt=testIBS,confInt=confInt,confLevel=confLevel)
    ## }
    # }}}
    # {{{ looping output
    ##     if (multiSplitTest==TRUE)
    ##       loopOut=list(PredCindexStepB=PredCindexStepB,testedResid=testedResid)
    ##     else
    loopOut=list(PredCindexStepB=PredCindexStepB)
    ##     if (keepResiduals==TRUE)  
    ##       loopOut=c(loopOut,list(Residuals=lapply(Residuals,function(R){
    ##         R[,sindex(eval.times=testTimes,jump.times=times)]
    ##       })))

    if (!is.null(getFromModel)){
      loopOut=c(loopOut,list(ModelParameters=ModelParameters))
    }
    loopOut
  }
  ## })
  b <- 1
  ## if (require(foreach)){
  if (missing(slaveseed)||is.null(slaveseed))
    slaveseed <- sample(1:1000000,size=B,replace=FALSE)
  Looping <- foreach (b= 1:B) %dopar% step(b,slaveseed[[b]])
  ## }
  ## else{
  ## Looping <- lapply(1:B,function(b){step(b,seed=NULL)})
  ## }
  # }}}
  # {{{ output
  ## 
  ## 
  ##    1. a list of NF matrices each with B (rows) and NT columns
  ##       the prediction error curves
  ## 
  ## if (verbose==TRUE && B>1) cat("\n")
  BootstrapCrossValCindexMat <- lapply(1:NF,function(f){
    ## matrix with NT columns and b rows
    do.call("rbind",lapply(Looping,function(b){
      c.b <- b$PredCindexStepB[[f]]$Cindex.b
      c.b
      ## pairs.b <- b$PredCindexStepB[[f]]$Pairs.b
      ## conc.b <- b$PredCindexStepB[[f]]$Concordant.b
    }))
  })
  ## 
  ##    2. a list of NF average out-of-bag prediction error curves
  ##       with length NT
  ##
  BootstrapCrossValCindex <- lapply(BootstrapCrossValCindexMat,colMeans)
  out <- list(BootstrapCrossValCindex=BootstrapCrossValCindex)
  ## 
  ##   3. the results of B residual tests 
  ##
  ##   print(str(Looping))
  if (multiSplitTest==TRUE){
    out$testedResid <- lapply(Looping,function(x)x$testedResid)
  }
  ## 
  ##   4. model parameters
  ##
  if (!is.null(getFromModel)){
    out$ModelParameters <- lapply(1:NF,function(f){
      lapply(Looping,function(x)x$ModelParameters[[f]])
    })
  }
  ## 
  ##   5. bootstrap crossvalidation results
  ##
  if (keepMatrix==TRUE)
    out$BootstrapCrossValCindexMat <- BootstrapCrossValCindexMat
  ## 
  ##   6. residuals
  ##
  ##   if (keepResiduals==TRUE){
  ##     out$Residuals <- lapply(1:NF,function(f){
  ##       bootResiduals <- lapply(Looping,function(b){
  ##         b$Residuals[[f]]
  ##       })
  ##       names(bootResiduals) <- paste("testSample",1:B,sep=".")
  ##       bootResiduals
  ##     })
  ##     names(out$Residuals) <- names(object)
  ##   }
  out
  # }}}
}
