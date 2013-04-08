#########################################
# Function 'predictEventProb.selectFGR' #
#########################################

#Author: Rob C.M. van Kruijsdijk
#Date original version: 24-02-2013
#Author: Thomas A. Gerds
#Date current version: 06-04-2013
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
  crrstep.form <- as.formula(update(formula,paste(timevar,"~.")),env=NULL)
  capture.output(crrstep.fit <- do.call("crrstep",list(formula=crrstep.form,data=data,etype=Event,failcode=cause,cencode=cens.code,trace = FALSE,...)))
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

predictEventProb.selectFGR <- function(object,newdata,times,...){
  predictEventProb(object[[1]],newdata=newdata,times=times,...)
}

