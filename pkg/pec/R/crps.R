crps <- function(object,
                 models,
                 what,
                 times,
                 ## weights=NULL,
                 start){
  stopifnot(class(object)[1] == "pec")
  
  ## find the prediction models
  if (missing(models)) models <- 1:length(object$models)
  else
    if (!is.numeric(models))
      models <- names(object$models)[match(models,names(object$models))]
  
  ## times
  object.times <- object$time
  if(missing(times)) times <- object$maxtime
  if (any(times>object$maxtime)) {
    warning(paste("You asked to integrate until times where prediction error curves are not defined.", object$maxtime))
    times <- times[times<=object$maxtime]
  }
  if (!(object$exact || length(object.times)>100))
    warning("Only ", length(time)," time point",ifelse(length(times)==1,"","s")," used")
  
  ##  time range
  if (missing(start)) start <- object$start
  ## list of values 
  if (missing(what) || is.null(what)){
    what <- grep(c("Err$"),names(object),val=TRUE)
  }
  ##   if (length(weights)>0){
  ##     stopifnot(length(weights)==length(object$time))
  ##   }
  ##   else weigths <- NULL
  out <- lapply(what,function(w){
    est <- object[[w]][models]
    ##     if (length(weights>0)){
    ##       print("Weighted sum")
    ##       xmat <- lapply(xmat,function(u){u*weights})
    ## }
    y <- sapply(times,function(t){
      intx <- sapply(est, function(y){
        Dint(x=object.times,
             y=y,
             range=c(start,t))
      })
    })
    if (!is.null(dim(y))){
      tnames <- paste("time=",round(times,1),sep="")
      tnames[times<1] <- paste("time=",signif(times[times<1],2),sep="")
      colnames(y) <- paste("IBS[",start,";",tnames,"]",sep="")
      y}
  })
  NW <- length(what)
  NT <- length(times)
  if (NW==1)
    out <- out[[1]]
  else
    names(out) <- what
  if (NT==1){
    if(NW>1){
      out <- do.call("cbind",out)
      colnames(out) <- what
    }
  }
  out
}
## the name ibs is more intuitive for integrated Brier score
## whereas continuous ranked probability score is less well known
ibs <- crps
