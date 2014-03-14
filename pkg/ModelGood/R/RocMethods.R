# {{{ Roc.list
#' Comparing prediction models with Receiver operating characteristics and
#' Brier scores
#' 
#' Evaluating the performance of risk prediction models with binary status
#' response variable (case/control or similar). Roc curves are based on a
#' single marker, which is the probability predicted case probability by a
#' model. The area under the curve and the Brier score is used to summarize the
#' performance.
#' 
#' Bootstrap-crossvalidation techniques are implemented to estimate the
#' generalization performance of the model(s), i.e. the performance which can
#' be expected in new subjects.
#' 
#' Missing data in the response or in the input matrix cause a failure.
#' 
#' For each prediction model there must be a \code{predictStatusProb} method:
#' for example, to assess a prediction model which evaluates to a
#' \code{myclass} object one defines a function called
#' \code{predictStatusProb.myclass} with arguments
#' \code{object,newdata,cutpoints,train.data,...}, like this
#' 
#' myFit=myModel(Y~X,data=dat)
#' 
#' class(myFit)="myclass"
#' 
#' \code{predictStatusProb.myclass <-
#' function(object,newdata,cutpoints,train.data,...){ predict(object,
#' data=newdata,method="probabilities") out }}
#' 
#' Such a function takes the object which was fitted with train.data and
#' derives a matrix with predicted event status probabilities for each subject
#' in newdata (rows) and each cutpoint (column) of the response variable that
#' defines an event status.
#' 
#' Currently implemented are \code{predictStatusProb} methods for the following
#' R-functions: \describe{ \item{}{\code{glm} (from \code{library(stats)}}
#' \item{}{\code{lrm} (from \code{library(Design)}} \item{}{\code{rpart} (from
#' \code{library(rpart)})} \item{}{\code{randomForest} from
#' \code{library(randomForest)}} }
#' 
#' @aliases Brier Brier.list Brier.glm Brier.lrm Brier.randomForest Brier.rpart
#' Roc Roc.list Roc.glm Roc.lrm Roc.randomForest Roc.rpart
#' @usage
#' \method{Roc}{list} (object, formula, data, splitMethod="noSplitMethod",
#' noinf.method=c("simulate"), simulate="reeval", B, M, breaks, cbRatio=1,
#' RocAverageMethod="vertical",
#' RocAverageGrid=switch(RocAverageMethod, "vertical"=seq(0,1,.01),
#' "horizontal"=seq(1,0,-.01)), model.args=NULL, model.parms=NULL,
#' keepModels=FALSE, keepSampleIndex=FALSE, keepCrossValRes=FALSE,
#' keepNoInfSimu, slaveseed, cores=1, na.accept=0, verbose=FALSE, ...)
#' 
#' @param object A named list of prediction models. Each entry is either an
#' R-object for which a \link{predictStatusProb} method exists (see details) or
#' a \code{call} that can be evaluated to such an R-object.  For
#' cross-validation (e.g. when \code{splitMethod} is "boot632plus") all the
#' models in this list must include their \code{call} in their value.
#' @param formula A formula. If missing, use the formula of the (first) model
#' in object.  The left hand side is used to find the status response variable
#' in \code{data}.
#' @param data A data frame to validate the prediction models If missing, use
#' the data of the (first) model in object.
#' @param splitMethod Method for estimating the generalization error.
#' 
#' \code{none}:Assess the models in the same data where they are fitted. Yields
#' the apparent or re-substitution performance. Often overestimates the
#' generalization performance.
#' 
#' \code{BootCV}: Bootstrap cross validation. The prediction models are trained
#' on \code{B} bootstrap samples, that are either drawn with or without
#' replacement from \code{data} of the size \code{M}.  The model performance is
#' estimated with the observations that are NOT in the bootstrap sample.
#' 
#' \code{Boot632}: Linear combination of the apparent performance and the
#' Bootcv performance using the constant weight .632 (see references).
#' 
#' \code{Boot632plus}: Linear combination of apparent performance and Bootcv
#' using weights dependent on how the models perform in permuted data (see
#' references).
#' 
#' \code{NoInf}: Assess the models in permuted data.
#' @param noinf.method Method for computing the no-information Roc curve.
#' @param simulate If equal to \code{"reeval"} then the models are re-evaluated
#' in the current permuted data for computing the no-information Roc curve.
#' @param B Number of bootstrap samples. The default depends on argument
#' \code{splitMethod}.  When \code{splitMethod in
#' c("Bootcv","Boot632","Boot632plus"} the default is 100. For
#' \code{splitMethod="cvK"} \code{B} is the number of cross-validation cycles,
#' and -- default is 1.  For \code{splitMethod="none"} \code{B} is the number
#' of bootstrap simulations e.g. to obtain bootstrap confidence limits --
#' default is 0.
#' @param M The size of the bootstrap samples for cross-validation without
#' replacement.
#' @param breaks Break points for computing the Roc curve. Defaults to
#' \code{seq(0,1,.01)} for the Roc.list method and to
#' \code{sort(unique(breaks))} for the default method.
#' @param cbRatio Experimental and not yet tested. Cost/benefit ratio. Default
#' value is to 1, meaning that misclassified cases are as bad as misclassified
#' controls.
#' @param RocAverageMethod Method for averaging ROC curves across data splits.
#' If \code{"horizontal"} average crossvalidated specificities for fixed sensitivity values,
#' specified in \code{RocAverageGrid}, otherwise, if  \code{"vertical"},
#' average crossvalidated specificities for fixed sensitivity values.
#' See Fawcett, T. (2006) for details.
#' @param RocAverageGrid Grid points for the averaging of Roc curves. A sequence of values at which to compute averages across the ROC curves obtained for different data splits during crossvalidation.
#' @param model.args List of extra arguments that can be passed to the
#' \code{predictStatusProb} methods. The list must have an entry for each entry
#' in \code{object}.
#' @param model.parms List with exactly one entry for each entry in
#' \code{object}.  Each entry names parts of the value of the fitted models
#' that should be extracted and added to the output (see value).
#' @param keepModels If \code{FALSE} keep only the names of the elements of
#' object.  If \code{"Call"} then keep the call of the elements of object.
#' Else, add the object as it is to the output.
#' @param keepSampleIndex Logical. If \code{FALSE} remove the cross-validation
#' index (which tells who was in the learn and who in the validation set) from
#' the output list which otherwise is included in the method part of the output
#' list.
#' @param keepCrossValRes Logical. If \code{TRUE} add all \code{B}
#' crossvalidation results to the output (see value). Defaults to \code{TRUE}.
#' @param keepNoInfSimu Logical. If \code{TRUE} add the \code{B} results in
#' permuted data (for no-information performance) to the output (see value).
#' Defaults to \code{FALSE}.
#' @param slaveseed Vector of seeds, as long as \code{B}, to be given to the
#' slaves in parallel computing to control the models build in crossvalidation loop.
#' @param cores Number of cores for parallel computing.
#' Passed as the value of the argument \code{mc.cores}
#' when calling \code{\link{mclapply}}.
#' @param na.accept Works only for "Bootcv" estimate of performance.  The
#' maximal number of failures during training the models to the bootstrap
#' samples. Usually a small number relative to \code{B}.
#' @param verbose if \code{TRUE} the procedure is reporting details of the
#' progress, e.g. it prints the current step in cross-validation procedures.
#' @param ... Difficult to explain
#' @return Object of class \code{Roc} or class \code{Brier} for which
#' \code{print}, \code{summary}, and \code{plot} methods are available.
#' 
#' The object includes the following components: \item{Roc}{ The Roc curve(s)
#' estimated according to the \code{splitMethod}. A matrix where each column
#' represents the estimated prediction error of a fit at the time points in
#' time.  } \item{AppRoc}{ The Roc curve(s) estimated when the model(s) are
#' evaluated in the same data where they were fitted. Only if
#' \code{splitMethod} is one of "NoInf", "Bootcv", "Boot632" or "Boot632plus",
#' since otherwise \code{repla} is "apparent" and then this is stored in
#' \code{Roc} as explained above.  } \item{BootcvRoc}{ The prediction error
#' when the model(s) are trained in the bootstrap sample and evaluated in the
#' data that are not in the bootstrap sample.  Only if \code{splitMethod} is
#' one of "Boot632" or "Boot632plus". When \code{splitMethod="Bootcv"} then the
#' \code{BootcvRoc} is stored in the component \code{PredRoc}.  }
#' \item{NoInfRoc}{ The prediction error when the model(s) are evaluated in the
#' permuted data.  Only if \code{splitMethod} is one of "Bootcv", "Boot632", or
#' "Boot632plus".  For \code{splitMethod="NoInf"} the \code{NoInfRoc} is stored
#' in the component \code{PredRoc}.  } \item{weight}{ The weight used to linear
#' combine the \code{AppRoc} and the \code{BootcvRoc} Only if
#' \code{splitMethod} is one of "Boot632", or "Boot632plus".  } \item{overfit}{
#' Estimated \code{overfit} of the model(s).  See references.  Only if
#' \code{splitMethod} is one of "Boot632", or "Boot632plus".  } \item{call}{The
#' call that produced the object} \item{models}{See keepModels}
#' \item{method}{The method used for estimation of the overfitting bias.}
#' @author Thomas Gerds \email{tag@@biostat.ku.dk}
#' @references Fawcett, T. (2006). An introduction to ROC analysis. Pattern
#' Recognition Letters, 27, 861-874.
#' 
#' Gerds, Cai & Schumacher (2008). The Performance of Risk Prediction Models.
#' Biometrical Journal, Vol 50, 4, 457-479.
#' 
#' Efron, Tibshirani (1997) Journal of the American Statistical Association 92,
#' 548--560 Improvement On Cross-Validation: The .632+ Bootstrap Method.
#' 
#' Wehberg, S and Schumacher, M (2004) A comparison of nonparametric error rate
#' estimation methods in classification problems. Biometrical Journal, Vol 46,
#' 35--47
#' @keywords models
##' @examples
##' 
##' ## Generate some data with binary response Y
##' ## depending on X1 and X2 and X1*X2
##' set.seed(40)
##' N <- 40
##' X1 <- rnorm(N)
##' X2 <- rbinom(N,1,.4)
##' expit <- function(x) exp(x)/(1+exp(x))
##' lp <- expit(1 + X1 + X2 - X1*X2)
##' Y <- factor(rbinom(N,1,lp))
##' dat <- data.frame(Y=Y,X1=X1,X2=X2)
##' 
##' ## fit a logistic model
##' lm1 <- glm(Y~X1,data=dat,family="binomial")
##' lm2 <- glm(Y~X1+X2,data=dat,family="binomial")
##' r1=Roc(list(lm1,lm2),cbRatio=1)
##' summary(r1)
##' Brier(list(lm1,lm2),cbRatio=1)
##' 
##' 
##' # crossing curves
##' set.seed(40)
##' N=40
##' Y=rbinom(N,1,.5)
##' X1=rnorm(N)
##' X1[Y==1]=rnorm(sum(Y==1),mean=rbinom(sum(Y==1),1,.5))
##' X2=rnorm(N)
##' X2[Y==0]=rnorm(sum(Y==0),mean=rbinom(sum(Y==0),1,.5))
##' dat <- data.frame(Y=Y,X1=X1,X2=X2)
##' lm1 <- glm(Y~X1,data=dat,family="binomial")
##' lm2 <- glm(Y~X2,data=dat,family="binomial")
##' plot(Roc(list(lm1,lm2),data=dat,verbose=0,cbRatio=1))
##' 
##' library(randomForest)
##' dat$Y=factor(dat$Y)
##' rf <- randomForest(Y~X2,data=dat)
##' rocCV=Roc(list(RandomForest=rf,LogisticRegression=lm2),
##'     data=dat,
##'     splitMethod="bootcv",
##'     B=3,
##'     cbRatio=1)
##'
##' bs <- Brier(list(LogisticRegression1=lm1,LogisticRegression=lm2),
##'     data=dat,
##'     splitMethod="boot632+",
##'     B=3,
##'     cbRatio=1)
##' plot(rocCV,diag=TRUE)
##' #'
#' @export
# {{{ UseMethod
Roc <- function(object,...){
    UseMethod("Roc",object=object)
}
# }}}
#' @S3method Roc list
Roc.list <- function(object,
                     formula,
                     data,
                     splitMethod="noSplitMethod",
                     noinf.method=c("simulate"),
                     simulate="reeval",
                     B,
                     M,
                     breaks,
                     cbRatio=1,
                     RocAverageMethod="vertical",
                     RocAverageGrid=switch(RocAverageMethod,"vertical"=seq(0,1,.01),"horizontal"=seq(1,0,-.01)),
                     model.args=NULL,
                     model.parms=NULL,
                     keepModels=FALSE,
                     keepSampleIndex=FALSE,
                     keepCrossValRes=FALSE,
                     keepNoInfSimu,
                     slaveseed,
                     cores=1,
                     na.accept=0,
                     verbose=FALSE,
                     ...){
    # }}}
    theCall=match.call()
    if (match("replan",names(theCall),nomatch=FALSE))
        stop("Argument name 'replan' has been replaced by 'splitMethod'.")
    # {{{ models
    NF <- length(object) 
    if (is.null(names(object)))names(object) <- sapply(object,function(o){
        class(o)[1]
    })
    else{names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])}
    object.names = names(object)
    names(object) <- make.names(names(object),unique=TRUE)
    # }}}
    # {{{ formula
    if (missing(formula)){
        if (class(object[[1]])=="formula")
            formula <- object[[1]]
        else
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
    if (SplitMethod$internal.name=="crossval")
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
        AppRoc <- Roc.default(object=pred,y=Y,breaks=breaks,cbRatio=cbRatio)
        AppAuc <- Auc.default(object=AppRoc$Sensitivity,Spec=AppRoc$Specificity)
        AppBS <- Brier.default(object=pred,y=Y,cbRatio=cbRatio)
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
                    innerNoInfRoc <- Roc.default(object=pred.index,y=data.index[,responseName],breaks=breaks,cbRatio=cbRatio)
                    innerNoInfBS <- Brier.default(object=pred.index,y=data.index[,responseName],cbRatio=cbRatio)
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
        # }}}
        # {{{ Bootcv aka BootstrapCrossValidation
        if (SplitMethod$internal.name %in% c("boot632plus","bootcv","boot632")){
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
                        innerBootcvRoc <- Roc.default(y=Y[vindex.index],pred.index,breaks=breaks,cbRatio=cbRatio)
                        innerBCVBS <- Brier.default(object=pred.index,y=Y[vindex.index],cbRatio=cbRatio)
                    }
                }
                list("innerBootcvRoc"=innerBootcvRoc,
                     "fit.parms"=fit.parms.index,
                     "failed"=failed,
                     "innerBCVBS"=innerBCVBS)
            }
            ## if (require(foreach)){
            b <- 1
            ## require(foreach)
            ## compute.BootcvRocList <- foreach::foreach (b = 1:B) %dopar% compute.step(runb=b,seed=slaveseed[[b]])
            compute.BootcvRocList <- parallel::mclapply(1:B,function(b){
                compute.step(runb=b,seed=slaveseed[[b]])
            },mc.cores=cores)
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
                  breaks=breaks,cbRatio=cbRatio))
    if (verbose) cat("\n")
    # }}}
    class(out) <- "Roc"
    out
}

# {{{ Roc.default, Roc.glm,etc.
#' @method Roc default
#' @S3method Roc default
Roc.default <- function(object,
                        y,
                        breaks,
                        cbRatio=1,
                        ## pv=FALSE,
                        ## confint=FALSE,
                        ## confint.method="exact",
                        keep.tables=FALSE,
                        keep.breaks=FALSE,...){
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
        count.TestPos <- N-prodlim::sindex(jump.times=object,eval.times=eval.times)
        count.TestNeg <- N-count.TestPos
        a <- count.DiseasePos-prodlim::sindex(jump.times=object[Disease==1],eval.times=eval.times)
        b <- count.DiseaseNeg-prodlim::sindex(jump.times=object[Disease==0],eval.times=eval.times)
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
    sens <- cbRatio*a/count.DiseasePos
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

#' @method Roc formula
#' @S3method Roc formula
Roc.formula <- function(object,formula,data,...){
    Roc.list(object=list("logistic.regression"=object),formula=object,data,...)
}

#' @method Roc glm
#' @S3method Roc glm
Roc.glm <- function(object,formula,data,...){
    stopifnot(object$family$family=="binomial")
    Roc.list(object=list(object),formula,data,...)
}

#' @method Roc lrm
#' @S3method Roc lrm
Roc.lrm <- function(object,formula,data,...){
    Roc.list(object=list(object),formula,data,...)
}

#' @method Roc ElasticNet
#' @S3method Roc ElasticNet
Roc.ElasticNet <- function(object,formula,data,...){
  Roc.list(list(object),formula,data,...)
}

#' @method Roc rpart
#' @S3method Roc rpart
Roc.rpart <- function(object,formula,data,...){
    Roc.list(object=list(object),formula,data,...)
}

#' @method Roc randomForest
#' @S3method Roc randomForest
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

