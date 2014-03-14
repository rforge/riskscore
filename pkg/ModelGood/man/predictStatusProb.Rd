\name{predictStatusProb}
\alias{predictStatusProb}
\alias{predictStatusProb.default}
\alias{predictStatusProb.glm}
\alias{predictStatusProb.lrm}
\alias{predictStatusProb.randomForest}
\alias{predictStatusProb.rpart}
\title{Probability Predictions}
\usage{
  \method{predictStatusProb}{glm}(object,newdata,...)
}
\arguments{
  \item{object}{A model from which predicted probabilities
  can be extracted for the indiviuals in newdata.}

  \item{newdata}{A data frame with new data for which the
  \code{object} should provide predict probabilities. In
  medical applications \code{newdata} will typically
  consist of the data of patients whose data were not used
  for building the model.}

  \item{...}{Additional arguments that are passed on to the
  current method.}
}
\value{
  A vector with the predicted status probability for each
  row in \code{NROW(newdata)}.
}
\description{
  This function can be used to extract probabilistic status
  predictions from various diagnostic and prognostic models
  with binary status response. The function invokes
  particular methods which depend on the 'class' of the
  first argument.
}
\details{
  The function delivers predicted probabilities tailored
  for the model performance measures of the package. These
  probabilities are extracted from a fitted model of class
  \code{CLASS} with the function
  \code{predictStatusProb.CLASS}.
}
\note{
  It is easy to write a new predictStatusProb method. For
  examples, see the existing methods (e.g. via
  \code{getAnywhere(predictStatusProb.glm)}.

  The performance, in particular when doing
  cross-validation where the model is evaluated many times,
  can be improved by supressing in the call to the model
  all the computations that are not needed for probability
  prediction, for example standard error calculations.
}
\author{
  Thomas A. Gerds \email{tag@biostat.ku.dk}
}
\seealso{
  \code{\link{predict}},\code{\link{Roc}}
}
\keyword{models}

