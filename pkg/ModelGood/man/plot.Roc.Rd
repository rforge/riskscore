\name{plot.Roc}
\alias{plot.Roc}
\title{ROC curves for risk prediction models}
\usage{
\method{plot}{Roc}(x, ylab = "Sensitivity", xlab = "1-Specificity", models,
  type = "l", shadow = FALSE, simu = FALSE, control, grid = FALSE,
  diag = TRUE, box = FALSE, lwd = 2, lty, col, add = FALSE,
  axes = TRUE, legend, auc, percent = TRUE, ...)
}
\arguments{
  \item{x}{object obtained with \code{Roc}}

  \item{ylab}{Label y-axis}

  \item{xlab}{Label x-axis}

  \item{models}{Selection of models to plot. Should be a
  subset of \code{names(x$models)}. Makes sense when
  \code{x} contains multiple ROC curves.}

  \item{type}{The line type}

  \item{shadow}{Experimental. Show results of
  cross-validation.}

  \item{simu}{Experimental. Show noinformation results.}

  \item{control}{Control which estimates of the ROC curves
  to draw.}

  \item{grid}{If \code{TRUE} add a grid in the background
  of the graph.}

  \item{diag}{If \code{TRUE} add a diagonal line.}

  \item{box}{If \code{TRUE} add a box around the graph.}

  \item{lwd}{Vector of line widths for the ROC curves.}

  \item{lty}{Vector of line types for the ROC curves.}

  \item{col}{Vector of colours for the ROC curves.}

  \item{add}{If \code{TRUE} add ROC curves to existing
  plot.}

  \item{axes}{If \code{TRUE} draw axes.}

  \item{legend}{If \code{TRUE} draw a legend.}

  \item{auc}{If \code{TRUE} add the area under the curve to
  the legend.}

  \item{percent}{If \code{TRUE} show percent axes.}

  \item{...}{Use for smart control of some plot elements.}
}
\value{
ROC curves
}
\description{
ROC curves for risk prediction models
}
\details{
Multiple ROC curves are shown in one graph.
}
\examples{
# generate som data
set.seed(40)
N=40
Y=rbinom(N,1,.5)
X1=rnorm(N)
X1[Y==1]=rnorm(sum(Y==1),mean=rbinom(sum(Y==1),1,.5))
X2=rnorm(N)
X2[Y==0]=rnorm(sum(Y==0),mean=rbinom(sum(Y==0),1,.5))
dat <- data.frame(Y=Y,X1=X1,X2=X2)

# fit two logistic regression models
lm1 <- glm(Y~X1,data=dat,family="binomial")
lm2 <- glm(Y~X2+X1,data=dat,family="binomial")
plot(Roc(list(lm1,lm2),data=dat))

# add the area under the curves

plot(Roc(list(lm1,lm2),data=dat),auc=TRUE)

# alternatively, one can directly work with formula objects:
plot(Roc(list(LR.X1=Y~X1,LR.X1.X2=Y~X2+X1),data=dat),auc=TRUE)

# beyond the logistic regression model.
# the following example is optimized for speed
# illustrating the syntax,
# and not for optimized for performance of the
# randomForest or elastic net
library(randomForest)
library(glmnet)
dat$Y=factor(dat$Y)
rf <- randomForest(Y~X1+X2,data=dat,ntree=10)
en <- ElasticNet(Y~X1+X2,data=dat,nfolds=10,alpha=0.1)
set.seed(6)
rocCV=Roc(list(RandomForest=rf,ElasticNet=en,LogisticRegression=lm2),
  data=dat,
  verbose=FALSE,
  splitMethod="bootcv",
  B=4,
  cbRatio=1)
plot(rocCV,yaxis.las=2,legend.title="4 bootstrap-crossvalidation steps")
}
\author{
Thomas A. Gerds <tag@biostat.ku.dk>
}
\seealso{
Roc
}

