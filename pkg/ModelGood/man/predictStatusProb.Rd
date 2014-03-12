\name{predictStatusProb}
\alias{predictStatusProb}
\alias{predictStatusProb.randomForest}
\alias{predictStatusProb.lrm}
\alias{predictStatusProb.default}
\alias{predictStatusProb.glm}
\alias{predictStatusProb.rpart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Probability Predictions}
\description{
  This function can be used to extract probabilistic status predictions 
  from various diagnostic and prognostic models with binary status
  response. The function invokes
  particular methods which depend on the 'class' of the first argument.
}
\usage{
 predictStatusProb(object, ...)
 \method{predictStatusProb}{glm}(object,newdata,...)
 \method{predictStatusProb}{randomForest}(object,newdata,...)
 \method{predictStatusProb}{lrm}(object,newdata,...)
 \method{predictStatusProb}{rpart}(object,newdata,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A model from which predicted probabilities can be
    extracted for the indiviuals in newdata.}
  \item{newdata}{A data frame containing the predictor variable combinations for
    which predictions are desired. In medical applications
    \code{newdata} will often consist of the data of patients that were
    not used for building the model.}
  \item{\dots}{Additional arguments that are passed on to the current method.
  }
}
\details{
  The function delivers predicted probabilities tailored for the model performance
  measures of the package. These probabilities are extracted from
  a fitted model of class \code{CLASS} with the function
  \code{predictStatusProb.CLASS}.
}

\value{
  A vector with the predicted status probability for each row in
  \code{NROW(newdata)}.
}
% \references{}

\author{Thomas A. Gerds \email{tag@biostat.ku.dk}}
\note{
  It is easy to write a new
  predictStatusProb method. For examples, see the
  existing methods (e.g. via \code{getAnywhere(predictStatusProb.glm)}.

  The performance, in particular
  when doing cross-validation
  where the model is evaluated many times, can be improved by
  supressing in the call to the model all the computations that are not
  needed for probability prediction, for example standard error
  calculations.
}

\seealso{\code{\link{predict}},\code{\link{Roc}}}
% \examples{
%  }

\keyword{models}% at least one, from doc/KEYWORDS
