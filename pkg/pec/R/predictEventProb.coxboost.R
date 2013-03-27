coxboost <- function(formula,data,cv=TRUE,cause,...){
  call <- match.call(expand.dots=FALSE)
  formula.names <- try(all.names(formula),silent=TRUE)
  if (!(formula.names[2]=="Hist")) stop("The left hand side of formula look like this: Hist(time,event).")
  actual.terms <- terms(formula,data=data)
  formula <- eval(call$formula)
  response <- model.response(model.frame(formula,data))
  Time <- as.numeric(response[,"time"])
  if (attr(response,"model")=="competing.risks"){
    ## adapt the event variable
    Event <- rep(2,NROW(response))
    thisCause <- as.numeric(response[,"event"]==cause)
    Event[thisCause==1] <- 1
    Status <- as.numeric(response[,"status"])
    Event[Status==0] <- 0
  }
  else{
    ## survival
    Event <- as.numeric(response[,"status"])
  }
  X <- model.matrix(actual.terms,data=data)[,-c(1),drop=FALSE]## remove intercept
  if (NCOL(X)<=1) stop("CoxBoost needs at least two covariates.")
  cv.defaults=list(maxstepno=200,K=10)
  CoxBoost.defaults=list(stepno=100)
  args <- SmartControl(call= list(...),
                       keys=c("cv","CoxBoost"),
                       ignore=c("formula","data","cv","cause"),
                       forced=list("cv"=list(time=Time,status=Event,x=X),"CoxBoost"=list(time=Time,status=Event,x=X)),
                       defaults=list("cv"=cv.defaults,"CoxBoost"=CoxBoost.defaults),
                       ignore.case=FALSE,
                       replaceDefaults=FALSE,
                       verbose=TRUE)
  if (cv==TRUE){
    cv.step <- do.call("cv.CoxBoost",args$cv)
    args$CoxBoost$stepno <- cv.step$optimal.step
  }
  cb <- do.call("CoxBoost",args$CoxBoost)
  out <- list(coxboost=cb,
              stepno=args$CoxBoost$stepno,
              call=call,
              formula=formula,
              response=response)
  class(out) <- "coxboost"
  out
}

predictSurvProb.coxboost <- function(object,newdata,times,...) {
  newcova <- model.matrix(terms(object$formula,data=newdata),
                          data=model.frame(object$formula,data=newdata,na.action=na.fail))[,-c(1)]
  newcova <- newcova[,object$coxboost$xnames]
  p <- predict(object$coxboost,newcova,type="risk",times=times)
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}


predictEventProb.coxboost <- function(object,newdata,times,cause,...){
  if (missing(cause)) stop("missing cause")
  if (attr(object$response,"model")!="competing.risks") stop("Not a competing risk object")
  newcova <- model.matrix(terms(object$formula,data=newdata),
                          data=model.frame(object$formula,data=newdata,na.action=na.fail))[,-c(1)]
  newcova <- newcova[,object$coxboost$xnames]
  p <- predict(object$coxboost,newdata=newcova,type="CIF",times=times)
  if (is.null(dim(p))) {
    if (length(p)!=length(times))
      stop("Prediction failed (wrong number of times)")
  }
  else{
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop("Prediction failed")
  }
  p
}
predictLifeYearsLost.rfsrc <- function(object, newdata, times, cause, ...){
  if (missing(cause)) stop("missing cause")
  cif <- predict(object,newdata=newdata,importance="none",...)$cif[,,cause,drop=TRUE]
  pos <- sindex(jump.times=object$time.interest,eval.times=times)
  lyl <- matrix(unlist(lapply(1:length(pos), function(j) {
    pos.j <- 1:(pos[j]+1)
    p <- cbind(0,cif)[,pos.j,drop=FALSE]
    time.diff <- diff(c(0, object$time.interest)[pos.j])
    apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
  })), ncol = length(pos))
  if (NROW(lyl) != NROW(newdata) || NCOL(lyl) != length(times))
    stop("Prediction of life-years-lost failed")
  lyl
}