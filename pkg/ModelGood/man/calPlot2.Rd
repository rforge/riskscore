\name{calPlot2}
\alias{calPlot2}
\title{Calibration plots for binary data}
\usage{
calPlot2(object, formula, data, splitMethod = "none", B = 1, M, showY,
  method = "nne", round = TRUE, bandwidth = NULL, q = 10,
  density = 55, add = FALSE, diag = !add, legend = !add, axes = !add,
  xlim, ylim, xlab = "Predicted event probability",
  ylab = "Observed proportion", col, lwd, lty, pch, cause = 1,
  percent = TRUE, giveToModel = NULL, na.action = na.fail, cores = 1,
  verbose = FALSE, ...)
}
\arguments{
  \item{object}{A named list of prediction models, where
  allowed entries are (1) R-objects for which a
  \link{predictStatusProb} method exists (see details), (2)
  a \code{call} that evaluates to such an R-object (see
  examples), (3) a matrix with predicted probabilities
  having as many rows as \code{data} in one column. For
  cross-validation all objects in this list must include
  their \code{call}.}

  \item{formula}{A survival or event history formula. The
  left hand side is used to compute the expected event
  status. If \code{formula} is \code{missing}, try to
  extract a formula from the first element in object.}

  \item{data}{A data frame in which to validate the
  prediction models and to fit the censoring model.  If
  \code{data} is missing, try to extract a data set from
  the first element in object.}

  \item{splitMethod}{Defines the internal validation
  design:

  \code{none/noPlan}: Assess the models in the give
  \code{data}, usually either in the same data where they
  are fitted, or in independent test data.

  \code{BootCv}: Bootstrap cross validation. The prediction
  models are trained on \code{B} bootstrap samples, that
  are either drawn with replacement of the same size as the
  original data or without replacement from \code{data} of
  the size \code{M}.  The models are assessed in the
  observations that are NOT in the bootstrap sample.}

  \item{B}{The number of cross-validation steps.}

  \item{M}{The size of the subsamples for
  cross-validation.}

  \item{showY}{If \code{TRUE} the observed data are shown
  as dots on the plot.}

  \item{method}{The method for estimating the calibration
  curve(s):

  \code{"nne"}: The expected event status is obtained in
  the nearest neighborhood around the predicted event
  probabilities.

  \code{"quantile"}: The expected event status is obtained
  in groups defined by quantiles of the predicted event
  probabilities.}

  \item{round}{If \code{TRUE} predicted probabilities are
  rounded to two digits before smoothing. This may have a
  considerable effect on computing efficiency in large data
  sets.}

  \item{bandwidth}{The bandwidth for \code{method="nne"}}

  \item{q}{The number of quantiles for
  \code{method="quantile"}.}

  \item{density}{Gray scale for observations.}

  \item{add}{If \code{TRUE} the line(s) are added to an
  existing plot.}

  \item{diag}{If \code{FALSE} no diagonal line is drawn.}

  \item{legend}{If \code{FALSE} no legend is drawn.}

  \item{axes}{If \code{FALSE} no axes are drawn.}

  \item{xlim}{Limits of x-axis.}

  \item{ylim}{Limits of y-axis.}

  \item{xlab}{Label for y-axis.}

  \item{ylab}{Label for x-axis.}

  \item{col}{Vector with colors, one for each element of
  object. Passed to \code{\link{lines}}.}

  \item{lwd}{Vector with line widths, one for each element
  of object. Passed to \code{\link{lines}}.}

  \item{lty}{lwd Vector with line style, one for each
  element of object. Passed to \code{\link{lines}}.}

  \item{pch}{Passed to \code{\link{points}}.}

  \item{cause}{For competing risks models, the cause of
  failure or event of interest}

  \item{percent}{If TRUE axes labels are multiplied by 100
  and thus interpretable on a percent scale.}

  \item{giveToModel}{List of with exactly one entry for
  each entry in \code{object}. Each entry names parts of
  the value of the fitted models that should be extracted
  and added to the value.}

  \item{na.action}{Passed to \code{\link{model.frame}}}

  \item{cores}{Number of cores for parallel computing.
  Passed as the value of the argument \code{mc.cores} when
  calling \code{\link{mclapply}}.}

  \item{verbose}{if \code{TRUE} report details of the
  progress, e.g. count the steps in cross-validation.}

  \item{...}{Used to control the subroutines: plot, axis,
  lines, legend. See \code{\link{SmartControl}}.}
}
\value{
list with elements: time, Frame and bandwidth (NULL for
method quantile).
}
\description{
Calibration plots for risk prediction models in for a
binary endpoint
}
\details{
For method "nne" the optimal bandwidth with respect to is
obtained with the function \code{\link{dpik}} from the
package \code{KernSmooth} for a box kernel function.
}
\examples{
set.seed(40)
N=40
Y=rbinom(N,1,.5)
X1=rnorm(N)
X1[Y==1]=rnorm(sum(Y==1),mean=rbinom(sum(Y==1),1,.5))
X2=rnorm(N)
X2[Y==0]=rnorm(sum(Y==0),mean=rbinom(sum(Y==0),3,.5))
dat <- data.frame(Y=Y,X1=X1,X2=X2)
lm1 <- glm(Y~X1,data=dat,family="binomial")
lm2 <- glm(Y~X2,data=dat,family="binomial")
calPlot2(list(lm1,lm2),data=dat)
}
\author{
Thomas Alexander Gerds
}
\references{
TA Gerds, PA Andersen, and Kattan MW. Calibration plots for
risk prediction models in the presence of competing risks.
Statistics in Medicine, page to appear, 2014.
}

