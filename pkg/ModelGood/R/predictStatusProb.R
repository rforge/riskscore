predictStatusProb <- function(object,...){
  UseMethod("predictStatusProb",object)
}

predictStatusProb.numeric <- function(object,newdata,...){
  stopifnot(NROW(object)==NROW(newdata))
  object
}


predictStatusProb.double <- function(object,newdata,...){
  stopifnot(NROW(object)==NROW(newdata))
  object
}

## predictStatusProb.NULL <- function(object,newdata,...){
## runif(NROW(newdata))
## }

predictStatusProb.glm <- function(object,newdata,...){
  if (object$family$family=="binomial")
    p <- predict(object,newdata=newdata,type="response")
  else{ stop("Currently only the binomial family is implemented for predicting a status from a glm object.")
      }
  p
}

predictStatusProb.BinaryTree <- function(object,newdata,...){
  treeresponse <- party::treeresponse
  P <- sapply(treeresponse(object,newdata=newdata),function(x)x[1])
  P
}

predictStatusProb.lrm <- function(object,newdata,...){
  P <- predict(object,newdata=newdata,type="fitted")
  P
}


predictStatusProb.rpart <- function(object,newdata,...){
  P <- as.numeric(predict(object,newdata=newdata,type="prob")[,2,drop=TRUE])
  P
}


predictStatusProb.randomForest <- function(object,newdata,...){
  stopifnot(!missing(newdata))
  P <- as.numeric(predict(object,newdata=newdata,type="prob")[,2,drop=TRUE])
  P
}

predictStatusProb.rfsrc <- function(object, newdata, times, ...){
  p <- predict(object,newdata=newdata,importance="none",...)$predicted[,2]
  p
}



