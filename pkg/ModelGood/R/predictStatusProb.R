#' Probability Predictions
#' 
#' This function can be used to extract probabilistic status predictions from
#' various diagnostic and prognostic models with binary status response. The
#' function invokes particular methods which depend on the 'class' of the first
#' argument.
#' 
#' The function delivers predicted probabilities tailored for the model
#' performance measures of the package. These probabilities are extracted from
#' a fitted model of class \code{CLASS} with the function
#' \code{predictStatusProb.CLASS}.
#' 
#' @aliases predictStatusProb predictStatusProb.randomForest
#' predictStatusProb.lrm predictStatusProb.default predictStatusProb.glm
#' predictStatusProb.rpart
#' @usage
#' \method{predictStatusProb}{glm}(object,newdata,...)
#' @param object A model from which predicted probabilities can be
#' extracted for the indiviuals in newdata.
#' @param newdata A data frame with new data for which the \code{object} should provide predict probabilities. In medical
#' applications \code{newdata} will typically consist of the data of 
#' patients whose data were not used for building the model.
#' @param ... Additional arguments that are passed on to the current
#' method.
#' @return A vector with the predicted status probability for each row in
#' \code{NROW(newdata)}.
#' @note It is easy to write a new predictStatusProb method. For examples, see
#' the existing methods (e.g. via \code{getAnywhere(predictStatusProb.glm)}.
#' 
#' The performance, in particular when doing cross-validation where the model
#' is evaluated many times, can be improved by supressing in the call to the
#' model all the computations that are not needed for probability prediction,
#' for example standard error calculations.
#' @usage predictStatusProb(object,newdata, ...)
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{predict}},\code{\link{Roc}}
#' @keywords models
#' @export predictStatusProb
predictStatusProb <- function(object,...){
  UseMethod("predictStatusProb")
}

##' @S3method predictStatusProb numeric
predictStatusProb.numeric <- function(object,newdata,...){
  stopifnot(NROW(object)==NROW(newdata))
  object
}

##' @S3method predictStatusProb formula
predictStatusProb.formula <- function(object,newdata,...){
  fit <- glm(object,data=newdata,family="binomial")
  predictStatusProb(fit,newdata=newdata,...)
}

##' @S3method predictStatusProb double
predictStatusProb.double <- function(object,newdata,...){
  stopifnot(NROW(object)==NROW(newdata))
  object
}

## predictStatusProb.NULL <- function(object,newdata,...){
## runif(NROW(newdata))
## }
##' @S3method predictStatusProb glm
predictStatusProb.glm <- function(object,newdata,...){
  if (object$family$family=="binomial")
    p <- predict(object,newdata=newdata,type="response")
  else{ stop("Currently only the binomial family is implemented for predicting a status from a glm object.")
      }
  p
}
##' @S3method predictStatusProb BinaryTree
predictStatusProb.BinaryTree <- function(object,newdata,...){
  treeresponse <- party::treeresponse
  P <- sapply(treeresponse(object,newdata=newdata),function(x)x[1])
  P
}

##' @S3method predictStatusProb lrm
predictStatusProb.lrm <- function(object,newdata,...){
  P <- predict(object,newdata=newdata,type="fitted")
  P
}

##' @S3method predictStatusProb rpart
predictStatusProb.rpart <- function(object,newdata,...){
  P <- as.numeric(predict(object,newdata=newdata,type="prob")[,2,drop=TRUE])
  P
}

##' @S3method predictStatusProb randomForest
predictStatusProb.randomForest <- function(object,newdata,...){
  stopifnot(!missing(newdata))
  P <- as.numeric(predict(object,newdata=newdata,type="prob")[,2,drop=TRUE])
  P
}

##' @S3method predictStatusProb rfsrc
predictStatusProb.rfsrc <- function(object, newdata, times, ...){
  p <- predict(object,newdata=newdata,importance="none",...)$predicted[,2]
  p
}



