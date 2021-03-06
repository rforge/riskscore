#' Calibration plots for right censored data
#' 
#' Calibration plots for risk prediction models in right censored survival and
#' competing risks data
#' 
#' For method "nne" the optimal bandwidth with respect to is obtained with the
#' function \code{\link{dpik}} from the package \code{KernSmooth} for a box
#' kernel function.
#' 
#' @param object A named list of prediction models, where allowed
#' entries are (1) R-objects for which a \link{predictSurvProb} method
#' exists (see details), (2) a \code{call} that evaluates to such an
#' R-object (see examples), (3) a matrix with predicted probabilities
#' having as many rows as \code{data} and as many columns as
#' \code{times}. For cross-validation all objects in this list must
#' include their \code{call}.
#' @param time The evaluation time point at predicted event
#' probabilities are plotted against pseudo-observed event status.
#' @param formula A survival or event history formula. The left hand
#' side is used to compute the expected event status. If
#' \code{formula} is \code{missing}, try to extract a formula from the
#' first element in object.
#' @param data A data frame in which to validate the prediction models
#' and to fit the censoring model. If \code{data} is missing, try to
#' extract a data set from the first element in object.
#' @param splitMethod Defines the internal validation design:
#' 
#' \code{none/noPlan}: Assess the models in the give \code{data}, usually
#' either in the same data where they are fitted, or in independent test data.
#' 
#' \code{BootCv}: Bootstrap cross validation. The prediction models
#' are trained on \code{B} bootstrap samples, that are either drawn
#' with replacement of the same size as the original data or without
#' replacement from \code{data} of the size \code{M}.  The models are
#' assessed in the observations that are NOT in the bootstrap sample.
#' @param B The number of cross-validation steps.
#' @param M The size of the subsamples for cross-validation.
#' @param outcome The method for estimating expected event status:
#' 
#' \code{"pseudo"}: Use average pseudo-values.  \code{"prodlim"}: Use
#' the product-limit estimate, i.e., apply the Kaplan-Meier method for
#' right censored survival and the Aalen-Johansen method for right
#' censored competing risks data.
#' @param showPseudo If \code{TRUE} and \code{outcome=="pseudo"} the
#' pseudo-values are shown as dots on the plot.
#' @param pseudo.col Colour for pseudo-values.
#' @param pseudo.pch Dot type (see par) for pseudo-values.
#' @param method The method for estimating the calibration curve(s):
#' 
#' \code{"nne"}: The expected event status is obtained in the nearest
#' neighborhood around the predicted event probabilities.
#' 
#' \code{"quantile"}: The expected event status is obtained in groups
#' defined by quantiles of the predicted event probabilities.
#' @param round If \code{TRUE} predicted probabilities are rounded to
#' two digits before smoothing. This may have a considerable effect on
#' computing efficiency in large data sets.
#' @param bandwidth The bandwidth for \code{method="nne"}
#' @param q The number of quantiles for \code{method="quantile"}.
#' @param jack.density Gray scale for pseudo-observations.
#' @param add If \code{TRUE} the line(s) are added to an existing
#' plot.
#' @param diag If \code{FALSE} no diagonal line is drawn.
#' @param legend If \code{FALSE} no legend is drawn.
#' @param axes If \code{FALSE} no axes are drawn.
#' @param xlim Limits of x-axis.
#' @param ylim Limits of y-axis.
#' @param xlab Label for y-axis.
#' @param ylab Label for x-axis.
#' @param col Vector with colors, one for each element of
#' object. Passed to \code{\link{lines}}.
#' @param lwd Vector with line widths, one for each element of
#' object. Passed to \code{\link{lines}}.
#' @param lty lwd Vector with line style, one for each element of
#' object.  Passed to \code{\link{lines}}.
#' @param pch Passed to \code{\link{points}}.
#' @param cause For competing risks models, the cause of failure or
#' event of interest
#' @param percent If TRUE axes labels are multiplied by 100 and thus
#' interpretable on a percent scale.
#' @param giveToModel List of with exactly one entry for each entry in
#' \code{object}. Each entry names parts of the value of the fitted
#' models that should be extracted and added to the value.
#' @param na.action Passed to \code{\link{model.frame}}
#' @param cores Number of cores for parallel computing.  Passed as
#' value of argument \code{mc.cores} to \code{\link{mclapply}}.
#' @param verbose if \code{TRUE} report details of the progress,
#' e.g. count the steps in cross-validation.
#' @param ... Used to control the subroutines: plot, axis, lines,
#' legend. See \code{\link{SmartControl}}.
#' @return list with elements: time, pseudoFrame and bandwidth (NULL for method
#' quantile).
#' @keywords survival 
##' @examples
##' 
##' library(prodlim)
##' library(lava)
##' library(riskRegression)
##' library(survival)
##' set.seed(13)
##' m <- crModel()
##' regression(m, from = "X1", to = "eventtime1") <- 1
##' regression(m, from = "X2", to = "eventtime1") <- 1
##' m <- addvar(m,c("X3","X4","X5"))
##' distribution(m, "X1") <- binomial.lvm()
##' distribution(m, "X4") <- binomial.lvm()
##' d1 <- sim(m,100)
##' d2 <- sim(m,100)
##' csc <- CSC(Hist(time,event)~X1+X2+X3+X4+X5,data=d1)
##' fgr <- FGR(Hist(time,event)~X1+X2+X3+X4+X5,data=d1,cause=1)
##' predict.crr <- cmprsk:::predict.crr
##' par(mar=c(5,5,5,5),cex=1.3)
##' calPlot(list(csc,fgr),
##'         time=5,
##'         legend.x=-0.3,
##'         legend.y=1.35,
##'         ylab="Observed event status",
##'         legend.legend=c("Cause-specific Cox regression","Fine-Gray regression"),
##'         legend.xpd=NA)
#' @author Thomas Alexander Gerds \email{tag@@biostat.ku.dk}
#' @export calPlot
calPlot <- function(object,
                    time,
                    formula,
                    data,
                    splitMethod="none",
                    B=1,
                    M,
                    outcome=c("pseudo","prodlim"),
                    showPseudo,
                    pseudo.col=NULL,
                    pseudo.pch=NULL,
                    method="nne",
                    round=TRUE,
                    bandwidth=NULL,
                    q=10,
                    jack.density=55,
                    add=FALSE,
                    diag=!add,
                    legend=!add,
                    axes=!add,
                    xlim=c(0,1),
                    ylim=c(0,1),
                    xlab = "Predicted event probability",
                    ylab = "Pseudo-observed event status",
                    col,
                    lwd,
                    lty,
                    pch,
                    cause=1,
                    percent=TRUE,
                    giveToModel=NULL,
                    na.action=na.fail,
                    cores=1,
                    verbose=FALSE,
                    ...){

    if (missing(showPseudo))
        showPseudo <- ifelse(add||(outcome!="pseudo"),FALSE,TRUE)
    
    # {{{ find number of objects and lines
    cobj=class(object)[[1]]
    if (cobj!="list"){
        object <- list(object)
    }
    if (is.null(names(object))){
        names(object) <- sapply(object,function(o)class(o)[1])
        names(object) <- make.names(names(object),unique=TRUE)
    }
    else{
        names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
    }
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
        if (length(grep("~",as.character(object[[1]]$call$formula)))==0){
            stop(paste("Argument formula is missing and first model has no usable formula:",as.character(object[[1]]$call$formula)))
        } else{
            ftry <- try(formula <- eval(object[[1]]$call$formula),silent=TRUE)
            if ((class(ftry)=="try-error") || match("formula",class(formula),nomatch=0)==0)
                stop("Argument formula is missing and first model has no usable formula.")
            else if (verbose)
                warning("Formula missing. Using formula from first model")
        }
    }
    
    m <- model.frame(formula,data,na.action=na.action)
    response <- model.response(m)
    if (match("Surv",class(response),nomatch=FALSE))
        model.type <- "survival"
    else
        model.type <- attr(response,"model")
    if (is.null(model.type) & length(unique(response))==2)
        model.type <- "binary"
    if (!model.type=="binary"){
        neworder <- order(response[,"time"],-response[,"status"])
        response <- response[neworder,,drop=FALSE]
        Y <- response[,"time"]
        status <- response[,"status"]
        data <- data[neworder,]
        # }}}
        # {{{ prediction timepoint 

        if (missing(time))
            time <- median(Y)
        else
            if (length(time)>1)
                stop("Please specify only one time point.")
    }

    # }}}
    # {{{ compute pseudo-values

    #  require(pseudo)
    #  jack=pseudosurv(time=Y,event=status,tmax=time)[[3]]
    predictHandlerFun <- switch(model.type,
                                "binary"="predictStatusProb",
                                "competing.risks"="predictEventProb",
                                "survival"="predictSurvProb")
    if (model.type=="binary")
        if (is.factor(response))
            jack <- as.numeric(response==levels(response)[2])
        else
            jack <- as.numeric(response)
    ## ==levels(response)[1])
    else{
        margForm <- reformulate("1",response=formula[[2]])
        margFit <- prodlim::prodlim(margForm,data=data)
        jack <- prodlim::jackknife(margFit,cause=cause,times=time)
    }

    # }}}
    # {{{ call smartControls
    axis1.DefaultArgs <- list(side=1,las=1,at=seq(0,ylim[2],ylim[2]/4))
    ## if (showPseudo==TRUE){
    ## at2 <- seq(0,1,.25)
    ## if (min(jack)<0) at2 <- c(round(min(jack),2),at2)
    ## if (max(jack)>1) at2 <- c(at2,round(max(jack),2))
    ## axis2.DefaultArgs <- list(side=2,las=2,at=at2,mgp=c(4,1,0))
    ## }
    ## else{
    axis2.DefaultArgs <- list(side=2,las=2,at=seq(0,ylim[2],ylim[2]/4),mgp=c(4,1,0))
    ## }
    legend.DefaultArgs <- list(legend=names(object),lwd=lwd,col=col,lty=lty,cex=1.5,bty="n",y.intersp=1.3,x="topleft")
    lines.DefaultArgs <- list(type="l")
    if (missing(ylim)){
        if (showPseudo){
            ylim <- c(min(jack),max(jack))
        }
        else
            ylim <- c(0,1)
    }
    if (missing(xlim)){
        xlim <- c(0,1)
    }
    plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,ylab="",xlab=xlab)
    smartA <- prodlim::SmartControl(call= list(...),keys=c("plot","lines","legend","axis1","axis2"),ignore=NULL,ignore.case=TRUE,defaults=list("plot"=plot.DefaultArgs,"lines"=lines.DefaultArgs,"legend"=legend.DefaultArgs,"axis1"=axis1.DefaultArgs,"axis2"=axis2.DefaultArgs),forced=list("plot"=list(axes=FALSE),"axis1"=list(side=1)),verbose=TRUE)

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
            apppred <- switch(model.type,
                              "competing.risks"={do.call(predictHandlerFun,list(object[[f]],newdata=data,times=time,cause=cause))},
                              "survival"={do.call(predictHandlerFun,list(object[[f]],newdata=data,times=time))},
                              "binary"={do.call(predictHandlerFun,list(object[[f]],newdata=data))})
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
                switch(model.type,
                       "competing.risks"={do.call(predictHandlerFun,list(object=fit.k,newdata=val.k,times=time,cause=cause))},
                       "survival"={do.call(predictHandlerFun,c(list(object=fit.k,newdata=val.k,times=time),extraArgs))},
                       "binary"={do.call(predictHandlerFun,list(object=fit.k,newdata=val.k))})
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
        ## predframe <- na.omit(predframe)
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
        pred.list <- parallel::mclapply(1:B,function(b){
            if (verbose==TRUE) internalTalk(b,B)
            jackRefit <- FALSE
            vindex.b <- match(1:N,unique(ResampleIndex[,b]),nomatch=0)==0
            val.b <- data[vindex.b,,drop=FALSE]
            if (jackRefit){
                margFit.b <- prodlim::prodlim(margForm,data=val.b)
                jack.b <- prodlim::jackknife(margFit.b,cause=cause,times=time)
            }
            else{
                jack.b <- jack[match(1:N,unique(ResampleIndex[,b]),nomatch=0)==0]
            }
            train.b <- data[ResampleIndex[,b],,drop=FALSE]
            frame.b <- data.frame(jack=jack.b)
            bootpred <- do.call("cbind",lapply(1:NF,function(f){
                fit.b <- internalReevalFit(object=object[[f]],data=train.b,step=b,silent=FALSE,verbose=verbose)
                extraArgs <- giveToModel[[f]]
                try2predict <- try(pred.b <- switch(model.type,
                                                    "competing.risks"={do.call(predictHandlerFun,list(object=fit.b,newdata=val.b,times=time,cause=cause))},
                                                    "survival"={do.call(predictHandlerFun,c(list(object=fit.b,newdata=val.b,times=time),extraArgs))},
                                                    "binary"={do.call(predictHandlerFun,list(object=fit.b,newdata=val.b))}),silent=TRUE)
                if (inherits(try2predict,"try-error")==TRUE){
                    rep(NA,NROW(val.b))
                }else{
                    pred.b
                }
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
    method <- match.arg(method,c("quantile","nne"))
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
                       y <- unlist(predict(prodlim::prodlim(form.pcut,data=cbind(data,p=pcut)),
                                           cause=cause,
                                           newdata=data.frame(pcut=levels(pcut)),
                                           times=time,
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
                       if (round==TRUE){
                           if (!is.null(bandwidth) && bandwidth>=1){
                               message("No need to round predicted probabilities to calculate calibration in the large")
                           } else{
                               p <- round(p,2)
                           }
                       }
                       p <- na.omit(p)
                       if (no <- length(attr(p,"na.action")))
                           warning("calPlot: removed ",no," missing values in risk prediction.",call.=FALSE,immediate.=TRUE)
                       if (is.null(bandwidth)){
                           if (length(p)>length(apppred[,f+1])){
                               bw <- prodlim::neighborhood(apppred[,f+1])$bandwidth
                           }else{
                               bw <- prodlim::neighborhood(p)$bandwidth
                           }
                       } else{
                           bw <- bandwidth
                       }
                       if (bw>=1){
                           ## calibration in the large
                           plotFrame <- data.frame(x=mean(p),y=mean(jackF))
                       } else{
                           nbh <- prodlim::meanNeighbors(x=p,y=jackF,bandwidth=bw)
                           plotFrame <- data.frame(x=nbh$uniqueX,y=nbh$averageY)
                       }
                       attr(plotFrame,"bandwidth") <- bw
                       plotFrame
                   }else{
                       form.p <- reformulate("p",response=formula[[2]])
                       y <- unlist(predict(prodlim::prodlim(form.p,data=cbind(data,p=p)),
                                           cause=cause,
                                           newdata=data.frame(p=sort(p)),
                                           times=time,
                                           type=ifelse(model.type=="competing.risks","cuminc","surv")))
                       plotFrame <- data.frame(x=sort(p),y=y)
                   }
               })
    })
    # }}}
    # {{{ plot an empty frame
  
    if (add==FALSE){
        do.call("plot",smartA$plot)
        if (axes){
            if (percent){
                smartA$axis2$labels <- paste(100*smartA$axis2$at,"%")
                ## if (min(jack)<0) smartA$axis2$labels[1] <- ""
                ## if (max(jack)>1) smartA$axis2$labels[length(smartA$axis2$labels)] <- ""
                smartA$axis1$labels <- paste(100*smartA$axis1$at,"%")
            }
            do.call("axis",smartA$axis1)
            mgp2 <- smartA$axis2$mgp
            if (length(mgp2)>0){
                oldmgp <- par()$mgp
                par(mgp=mgp2)
                smartA$axis2 <- smartA$axis2[-match("mgp",names(smartA$axis2),nomatch=0)]
                title(ylab=ylab)
            }
            ## print(par()$mgp)
            do.call("axis",smartA$axis2)
            if (length(mgp2)>0){
                par(mgp=oldmgp)
            }
        }
    }
    if (diag){
        segments(x0=0,y0=0,x1=1,y1=1,col="gray77",lwd=2,xpd=FALSE)
    }
    ##   do.call("abline",c(list(a=0,b=1),list(col="gray77",lwd=2,xpd=FALSE)))
    # }}}
    # {{{ add lines and pseudovalues
    nix <- lapply(1:NF,function(f){
        if (is.null(pseudo.col)){
            ccrgb=as.list(col2rgb(col[f],alpha=TRUE))
            names(ccrgb) <- c("red","green","blue","alpha")
            ccrgb$alpha <- jack.density
            jack.col <- do.call("rgb",c(ccrgb,list(max=255)))
        }
        else
            jack.col <- pseudo.col
        if (is.null(pseudo.pch)) pseudo.pch <- 1
        if (showPseudo) {
            points(apppred[,f+1],apppred[,1],col=jack.col,pch=pseudo.pch)
        }
        plotFrame <- plotFrames[[f]]
        if(NROW(plotFrame)==1){
            plottype <- "p"
        } else{
            if (method=="quantile"){
                plottype <- "b"
            } else{
                plottype <- "l"
            }
        }
        with(na.omit(plotFrame),lines(x,y,col=col[f],lwd=lwd[f],lty=lty[f],type=plottype))
        ## }
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
    if (model.type=="binary")
        out <- list(pseudoFrame=plotFrames,bandwidth=sapply(plotFrames,function(x)attr(x,"bandwidth")))
    else
    out <- list(time=time,pseudoFrame=plotFrames,bandwidth=sapply(plotFrames,function(x)attr(x,"bandwidth")))
  invisible(out)
  # }}}
}
