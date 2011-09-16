
## rsf.default=randomSurvivalForest:::rsf.default 
predictSurvProb.rsf <- function (object, newdata, times, ...)  
{ 
  N <- NROW(newdata) 
  NT <- length(times)
  class(object) <- c("rsf", "grow")
  S <- exp(-predict.rsf(object, test=newdata)$ensemble)
  if (N==1) S <- matrix(S,nrow=1)
  Time <- object$timeInterest 
  p <- cbind(1,S)[,1+sindex(Time, times),drop=FALSE] 
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))  
    stop("Prediction failed") 
  p 
} 


#predictSurvProb.rsf <- function (object, newdata, times,...) {
#  require(randomSurvivalForest)
#  N <- NROW(newdata)
#  NT <- length(times)
#  class(object) <- c("rsf", "grow")
#  rsf.default <- randomSurvivalForest:::rsf.default
#  H <- predict.rsf(object=object, newdata=newdata)$ensemble
#  S <- exp(-H)
#  Time <- object$timeInterest
#  p <- cbind(1,S)[,1+sindex(Time, times), drop=FALSE]
#  print(dim(p))
#  print(c(NROW(newdata),length(times)))
#  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
#  stop("Prediction failed")
#  p
#}    

