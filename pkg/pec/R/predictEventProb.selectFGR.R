#########################################
# Function 'predictEventProb.selectFGR' #
#########################################

#Author: Rob C.M. van Kruijsdijk
#Date original version: 24-02-2013
#Contributor: Thomas A. Gerds
#Date current version: 06-04-2013


#' Stepwise variable selection in the Fine & Gray regression competing risk
#' model
#' 
#' This is a wrapper function which first selects variables in the Fine & Gray
#' regression model using \code{crrstep} from the \code{crrstep} package and
#' then returns a fitted Fine & Gray regression model with the selected
#' variables.
#' 
#' 
#' @param formula A formula whose left hand side is a \code{Hist} object -- see
#' \code{\link{Hist}}.  The right hand side specifies (a linear combination of)
#' the covariates. See examples below.
#' @param data A data.frame in which all the variables of \code{formula} can be
#' interpreted.
#' @param cause The failure type of interest. Defaults to \code{1}.
#' @param rule Rule to pass on to crrstep ("AIC", "BIC" or "BICcr"), also see
#' \code{crrstep}
#' @param direction see \code{crrstep}
#' @param \dots Further arguments passed to \code{crrstep}.
#' @author Rob C.M. van Kruijsdijk \email{R.C.M.vanKruijsdijk@@umcutrecht.nl}
#' 
#' Thomas Alexander Gerds \email{tag@@biostat.ku.dk}
#' @keywords survival
#' @examples
#' 
#' library(riskRegression)
#' library(prodlim)
#' library(pec)
#' library(cmprsk)
#' 
#' d <- SimCompRisk(100)
#' newvars <- as.data.frame(matrix(runif(8*100),nrow=100))
#' colnames(newvars) <- c('X3','X4','X5','X6','X7','X8','X9','X10')
#' d <- cbind(d,newvars)
#' require(cmprsk)
#' predict.crr <- cmprsk:::predict.crr
#' fg <- FGR(Hist(time, cause) ~ X1 + X2 + X3 + X4,cause=1,data=d)
#' sfg <- selectFGR(Hist(time, cause) ~ X1 + X2 + X3 + X4,
#' 		 cause=1,
#' 		 data=d,
#' 		 rule="AIC",
#' 		 direction="backward") 
#' predictEventProb(sfg,times=100,newdata=d)
#' 
#' ## apparent performance
#' pec(list(full.model=fg,selected.model=sfg),
#'     formula=Hist(time, cause)~1,
#'     data=d)
#' 
#' ## bootstrap cross-validation performance
#' set.seed(7)
#' pec(list(full.model=fg,selected.model=sfg),
#'     formula=Hist(time, cause)~1,
#'     data=d,
#'     B=5,
#'     splitMethod="bootcv")
#' 
#' @export selectFGR
selectFGR <- function(formula,data,cause=1,rule="AIC", direction="backward",...){
    if (!require(riskRegression)) stop("This function requires library riskRegression")
    if (!require(crrstep)) stop("This function requires library crrstep")
    if (missing(data)) stop("Argument 'data' is missing")
    if (missing(formula)) stop("Argument 'formula' is missing")
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    if (match("subset",names(call),nomatch=FALSE))
        stop("Subsetting of data is not possible.")
    m <- model.frame(formula,data)
    response <- model.response(m)
    cens.code <- as.numeric(attr(response,"cens.code"))
    timevar <- colnames(response)[1]
    if (attr(response,"model")=="competing.risks"){
        Event <- rep(2,NROW(response))
        thisCause <- as.numeric(response[,"event"]==cause)
        Event[thisCause==1] <- 1
        Status <- as.numeric(response[,"status"])
        Event[Status==0] <- 0
    }
    else{
        stop("This does not look like a competing risk setting.\nMaybe there is only one event type in the data?")
    }
    crrstep.form <- as.formula(update(formula,paste(timevar,"~1+.")))
    capture.output(crrstep.fit <- do.call("crrstep",list(formula=crrstep.form,
                                                         data=data,
                                                         etype=Event,
                                                         failcode=cause,
                                                         cencode=cens.code,
                                                         trace = FALSE,
                                                         ...)))
    if (length(crrstep.fit$coefficients)==0){
        newform <- as.formula(update(formula,.~1),env=NULL)
        newfit <- prodlim(newform,data=data)
    }
    else{
        newform <- as.formula(update(formula,paste(".~",paste(rownames(crrstep.fit$coefficients),collapse="+"))),env=NULL)
        newfit <- FGR(newform,data=data,cause=cause)
        newfit$call$formula <- newform
    }
    out <- list(fit=newfit,In=rownames(crrstep.fit$coefficients))
    out$call <- match.call()
    class(out) <- "selectFGR"
    out
}

##' @S3method predictEventProb selectFGR
predictEventProb.selectFGR <- function(object,newdata,times,...){
    predictEventProb(object[[1]],newdata=newdata,times=times,...)
}

