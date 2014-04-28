\name{print.pec}
\alias{print.pec}
\alias{summary.pec}
\title{Printing a `pec' (prediction error curve) object.}
\usage{
\method{print}{pec}(x, times, digits = 3, what = NULL, ...)
}
\arguments{
  \item{x}{Object of class \code{pec}}

  \item{object}{Object of class \code{pec}}

  \item{times}{Time points at which to show the values of
  the prediction error curve(s)}

  \item{what}{What estimate of the prediction error curve
  to show. Should be a string matching an element of the
  object. The default is determined by splitMethod.}

  \item{models}{Which models in the list
  \code{object$models} should be shown. Defaults to all
  models.}

  \item{digits}{Number of decimals used in tables.}

  \item{print}{Set to FALSE to suppress printing.}

  \item{\dots}{Not used}
}
\value{
The first argument in the invisible cloak.
}
\description{
Print the important arguments of call and the prediction
error values at selected time points.
}
\author{
Thomas A. Gerds \email{tag@biostat.ku.dk}
}
\seealso{
\code{\link{pec}}
}
\keyword{survival}

