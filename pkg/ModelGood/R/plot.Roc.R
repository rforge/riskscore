##' ROC curves for risk prediction models 
##'
##' Multiple ROC curves are shown in one graph.
##' @title ROC curves for risk prediction models 
##' @param x object obtained with \code{Roc}
##' @param ylab Label y-axis
##' @param xlab Label x-axis
##' @param models Selection of models to plot. Should be a subset of
##' \code{names(x$models)}. Makes sense when \code{x} contains
##' multiple ROC curves.
##' @param shadow Experimental. Show results of cross-validation.
##' @param simu Experimental. Show noinformation results.
##' @param control Control which estimates of the ROC curves to draw.
##' @param grid If \code{TRUE} add a grid in the background of the
##' graph.
##' @param diag If \code{TRUE} add a diagonal line.
##' @param box If \code{TRUE} add a box around the graph.
##' @param lwd Vector of line widths for the ROC curves.
##' @param lty Vector of line types for the ROC curves.
##' @param col Vector of colours for the ROC curves.
##' @param add If \code{TRUE} add ROC curves to existing plot.
##' @param axes If \code{TRUE} draw axes.
##' @param legend If \code{TRUE} draw a legend.
##' @param percent If \code{TRUE} show percent axes.
##' @param ... Use for smart control of some plot elements.
##' @return ROC curves
##' @seealso Roc
##' @examples
#' 
#' set.seed(40)
#' N=40
#' Y=rbinom(N,1,.5)
#' X1=rnorm(N)
#' X1[Y==1]=rnorm(sum(Y==1),mean=rbinom(sum(Y==1),1,.5))
#' X2=rnorm(N)
#' X2[Y==0]=rnorm(sum(Y==0),mean=rbinom(sum(Y==0),1,.5))
#' dat <- data.frame(Y=Y,X1=X1,X2=X2)
#' lm1 <- glm(Y~X1,data=dat,family="binomial")
#' lm2 <- glm(Y~X2,data=dat,family="binomial")
#' plot(Roc(list(lm1,lm2),data=dat,verbose=0,crRatio=1))
#' library(randomForest)
#' dat$Y=factor(dat$Y)
#' rf <- randomForest(Y~X2,data=dat)
#' set.seed(6)
#' rocCV=Roc(list(RandomForest=rf,LogisticRegression=lm2),
#'   data=dat,
#'   verbose=FALSE,
#'   splitMethod="bootcv",
#'   B=4,
#'   crRatio=1)
#' plot(rocCV,diag=TRUE,yaxis.las=2,legend.title="4 bootstrap-crossvalidation steps")
##' 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
##' @method plot Roc
##' @S3method plot Roc
plot.Roc <- function(x,
                     ylab="True detection rate",
                     xlab="False positive rate",
                     models,
                     shadow=FALSE,
                     simu=FALSE,
                     control,
                     grid=FALSE,
                     diag=FALSE,
                     box=FALSE,
                     lwd=2,
                     lty,
                     col,
                     add=FALSE,
                     axes=TRUE,
                     legend,
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
    do <- list("Roc"=1,"AppRoc"=0,"BootcvRoc"=0,"NoInfRoc"=0)
    ## if (roc==TRUE) do[["Roc"]] <- 1 
    ## if (apparent==TRUE) do[["AppRoc"]] <- do[["Roc"]] + 1
    ## if (BCV==TRUE) do[["BootcvRoc"]] <- do[["Roc"]] + (do[["AppRoc"]]>0) + 1
    ## if (noinf==TRUE) do[["NoInfRoc"]] <- do[["Roc"]] + (do[["AppRoc"]]>0) + (do[["BootcvRoc"]]>0) + 1
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
    xaxis.DefaultArgs <- list(at=seq(0,1,.25))
    yaxis.DefaultArgs <- list(at=seq(0,1,.25))
    smartA <- prodlim::SmartControl(call=  list(...),
                                    keys=c("plot","legend","xaxis","yaxis"),
                                    ignore=c("x","add","col","lty","lwd","ylim","xlim","xlab","ylab","legend","percent","grid","box","axes","roc","shadow","simu","BCV","apparent","noinf","control","diag"),
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
            if (percent & is.null(smartA$yaxis$labels))
                smartA$yaxis$labels <- paste(100*smartA$yaxis$at,"%")
            do.call("axis",smartA$yaxis)
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
                shadow.control <- c(shadow.control,list(type="l",col="gray77"))
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
                simu.control <- c(simu.control,list(type="l",col="gray77"))
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
                w.control <- c(w.control,list(type="l",lwd=lwd,lty=lty[[m]][w],col=col[[m]][w]))
                w.control <- w.control[!duplicated(names(w.control))]
                X <- 1-x[[W]][[m]]$Specificity
                Y <- x[[W]][[m]]$Sensitivity
                do.call("lines.default",c(list(X,Y), w.control))
            }
        }
    }
    # legend
    # --------------------------------------------------------------------
    if (missing(legend))
        if (length(names(x$models))>1)
            legend <- TRUE
        else
            legend <- FALSE
    if(legend==TRUE && !add && !is.null(names(x$models))){
        save.xpd <- par()$xpd
        par(xpd=TRUE)
        do.call("legend",smartA$legend)
        par(xpd=save.xpd)
    }
}

##' @S3method lines Roc
lines.Roc <- function(x,...){
    plot(x,add=TRUE,...)
}
