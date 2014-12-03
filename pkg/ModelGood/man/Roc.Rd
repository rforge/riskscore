\name{Roc}
\alias{Brier}
\alias{Brier.glm}
\alias{Brier.list}
\alias{Brier.lrm}
\alias{Brier.randomForest}
\alias{Brier.rpart}
\alias{Roc}
\alias{Roc.glm}
\alias{Roc.list}
\alias{Roc.lrm}
\alias{Roc.randomForest}
\alias{Roc.rpart}
\title{Comparing prediction models with Receiver operating characteristics and
Brier scores}
\usage{
\method{Roc}{list} (object, formula, data, splitMethod='noSplitMethod',
noinf.method=c('simulate'), simulate='reeval', B, M, breaks, cbRatio=1,
RocAverageMethod='vertical',
RocAverageGrid=switch(RocAverageMethod, 'vertical'=seq(0,1,.01),
'horizontal'=seq(1,0,-.01)), model.args=NULL, model.parms=NULL,
keepModels=FALSE, keepSampleIndex=FALSE, keepCrossValRes=FALSE,
keepNoInfSimu, slaveseed, cores=1, na.accept=0, verbose=FALSE, ...)
}
\arguments{
  \item{object}{A named list of R objects that represent
  predictive markers, prediction models, or prediction
  algorithms. The function \link{predictStatusProb} is
  called on the R objects to extract the predicted risk
  (see details). For cross-validation (e.g. when
  \code{splitMethod} is 'bootcv') all the R objects in this
  list must include a \code{call} which can be evaluated in
  a learning subset of the data.}

  \item{formula}{A formula whose left hand side is used to
  identify the binary outcome variable in \code{data}. If
  missing, use the formula of the (first) model in object.}

  \item{data}{A data set in which to validate the
  prediction models. If missing, the function tries to
  extract the data from the call of the (first) model in
  object.

  The data set needs to have the same structure, variable
  names, factor levels, etc., as the data in which the
  models were trained. If the subjects in data were not
  used to train the models given in \code{object}, this
  leads to an external validation situation.

  However, note that if one of the elements in
  \code{object} is a formula then it is evaluated in this
  data set.}

  \item{splitMethod}{Method for estimating the
  generalization error.

  \code{none}:Assess the models in the data given by
  \code{data}. If this data set coincides with the train
  data where the models were fitted this yields an apparent
  (or re-substitution) estimate of performance. Otherwise,
  this leads to an external validation situation.

  \code{bootCV}: Internal bootstrap cross validation. The
  prediction models are trained on \code{B} bootstrap
  samples of the \code{data}. Bootstrap samples are either
  drawn with replacement from \code{data} (same size), or
  without replacement of the size \code{M} where \code{M}
  is a number smaller than \code{NROW(data)}. The model
  performance parameters (Roc, Brier, AUC) are estimated
  with the observations that are NOT in the current
  bootstrap sample.

  \code{boot632}: Linear combination of the apparent
  performance and the BootCV performance using the constant
  weight .632 (see Efron & Tibshirani, 1997).

  \code{boot632plus}: Linear combination of apparent
  performance and Bootcv using weights dependent on how the
  models perform in permuted data (see Efron & Tibshirani,
  1997).

  \code{noinf}: Assess the models trained in permutations
  of \code{data}.}

  \item{noinf.method}{Experimental: For .632+ method the
  way to obtain no-information performance. This can either
  be 'simulate' or 'none'.}

  \item{simulate}{Experimental: For .632+ method. If
  \code{'reeval'} then the models are re-build in the
  current permuted data for computing the no-information
  Roc curve.}

  \item{B}{Number of repetitions for internal
  crossvalidation. The meaning depends on the argument
  \code{splitMethod}: When \code{splitMethod in
  c('Bootcv','Boot632','Boot632plus')} it is the number of
  bootstrap samples, default is 100. Otherwise it is
  ignored.}

  \item{M}{The size of the bootstrap samples for
  cross-validation without replacement.}

  \item{breaks}{Break points for computing the Roc curve.
  Defaults to \code{seq(0,1,.01)} when crossvalidation is
  applied, i.e., when \code{splitMethod in
  c('Bootcv','Boot632','Boot632plus')}. Otherwise use all
  unique values of the predictive marker.}

  \item{cbRatio}{Experimental. Cost/benefit ratio. Default
  value is 1, meaning that misclassified cases are as bad
  as misclassified controls.}

  \item{RocAverageMethod}{Method for averaging ROC curves
  across data splits. If \code{'horizontal'} average
  crossvalidated specificities for fixed sensitivity
  values, specified in \code{RocAverageGrid}, otherwise, if
  \code{'vertical'}, average crossvalidated specificities
  for fixed sensitivity values. See Fawcett, T. (2006) for
  details.}

  \item{RocAverageGrid}{Grid points for the averaging of
  Roc curves. A sequence of values at which to compute
  averages across the ROC curves obtained for different
  data splits during crossvalidation.}

  \item{model.args}{List of extra arguments that can be
  passed to the \code{predictStatusProb} methods. The list
  must have an entry for each entry in \code{object}.}

  \item{model.parms}{List with exactly one entry for each
  entry in \code{object}.  Each entry names parts of the
  value of the fitted models that should be extracted and
  added to the output (see value).}

  \item{keepModels}{If \code{FALSE} keep only the names of
  the elements of object.  If \code{'Call'} then keep the
  call of the elements of object. Else, add the object as
  it is to the output.}

  \item{keepSampleIndex}{Logical. If \code{FALSE} remove
  the cross-validation index (which tells who was in the
  learn and who in the validation set) from the output list
  which otherwise is included in the method part of the
  output list.}

  \item{keepCrossValRes}{Logical. If \code{TRUE} add all
  \code{B} crossvalidation results to the output (see
  value). Defaults to \code{TRUE}.}

  \item{keepNoInfSimu}{Logical. If \code{TRUE} add the
  \code{B} results in permuted data (for no-information
  performance) to the output (see value). Defaults to
  \code{FALSE}.}

  \item{slaveseed}{Vector of seeds, as long as \code{B}, to
  be given to the slaves in parallel computing to control
  the models build in crossvalidation loop.}

  \item{cores}{Number of cores for parallel computing.
  Passed as the value of the argument \code{mc.cores} when
  calling \code{\link{mclapply}}.}

  \item{na.accept}{For 'Bootcv' estimate of performance.
  The maximal number of bootstrap samples in which the
  training the models may fail This should usually be a
  small number relative to \code{B}.}

  \item{verbose}{if \code{TRUE} the procedure is reporting
  details of the progress, e.g. it prints the current step
  in cross-validation procedures.}

  \item{...}{Used to pass arguments to submodules.}
}
\value{
Object of class \code{Roc} or class \code{Brier}.

Depending on the \code{splitMethod} the object includes the
following components:

\item{Roc, Brier, Auc}{A list of Roc curve(s), Brier scores
(BS), and areas under the curves (Auc), one for each
element of argument \code{object}, estimated according to
\code{splitMethod}.}

\item{weight}{The weight used to linear combine the
\code{AppRoc} and the \code{BootcvRoc} Only available if
\code{splitMethod} is one of 'Boot632', or 'Boot632plus'. }

\item{overfit}{ Estimated \code{overfit} of the model(s).
Only if \code{splitMethod} is one of 'Boot632', or
'Boot632plus'.  }

\item{call}{The call that produced the object}

\item{models}{See keepModels}

\item{method}{Summary of the splitMethod used.}
}
\description{
Evaluation of the performance of risk prediction models
with binary status response variable (case/control or
similar). Roc curves are either based on a single
continuous marker, or on the probability prediction of an
event. Probability predictions are extracted from a given
(statistical) model, such as logistic regression, or
algorithm, such as random forest. The area under the curve
and the Brier score is used to summarize and compare the
performance.
}
\details{
All functions work on a list of models to ease comparison.

Bootstrap-crossvalidation techniques are implemented to
estimate the generalization performance of the model(s),
i.e., the performance which can be expected in new
subjects.

By default, when crossvalidation is involved, the ROC curve
is approximated on a grid of either sensitivities or
specificities and not computed at all unique changepoints
of the crossvalidated ROC curves, see Fawcett, T. (2006).
The (density of the) grid can be controlled with the
argument: RocAverageGrid

Missing data in the response or in the marker/predicted
risk cause a failure.

For each R object which potentially can predict a
probability for an event, there should be a corresponding
\code{predictStatusProb} method:

For example, to assess a prediction model which evaluates
to a \code{myclass} object one defines a function called
\code{predictStatusProb.myclass} with arguments
\code{object,newdata,...}. For example, the function
predictStatusProb.lrm looks like this:

predictStatusProb.lrm <- function(object,newdata,...){ p <-
as.numeric(predict(object,newdata=newdata,type='fitted'))
class(p) <- 'predictStatusProb' p }

Currently implemented are \code{predictStatusProb} methods
for the following R-functions: \describe{
\item{}{\code{numeric} (marker values are passed on)}
\item{}{\code{formula} (single predictor: extracted from
newdata and passed on, multiple predictors: projected to
score by logistic regression)} \item{}{\code{glm} (from
\code{library(stats)}} \item{}{\code{lrm} (from
\code{library(Design)}} \item{}{\code{rpart} (from
\code{library(rpart)})} \item{}{\code{BinaryTree} (from
\code{library(party)})} \item{}{\code{ElasticNet} (a
wrapper for glmnet from \code{library(glmnet)})}
\item{}{\code{randomForest} from
\code{library(randomForest)}} \item{}{\code{rfsrc} from
\code{library(randomForestSRC)}} }
}
\examples{
## Generate some data with binary response Y
## depending on X1 and X2 and X1*X2
set.seed(40)
N <- 40
X1 <- rnorm(N)
X2 <- abs(rnorm(N,4))
X3 <- rbinom(N,1,.4)
expit <- function(x) exp(x)/(1+exp(x))
lp <- expit(-2 + X1 + X2 + X3 - X3*X2)
Y <- factor(rbinom(N,1,lp))
dat <- data.frame(Y=Y,X1=X1,X2=X2)

# single markers, one by one
r1 <- Roc(Y~X1,data=dat)
plot(r1,col=1)
r2 <- Roc(Y~X2,data=dat)
lines(r2,col=2)

# or, directly multiple in one
r12 <- Roc(list(Y~X1,Y~X2),data=dat)
plot(r12)

## compare logistic regression
lm1 <- glm(Y~X1,data=dat,family="binomial")
lm2 <- glm(Y~X1+X2,data=dat,family="binomial")
r1=Roc(list(LR.X1=lm1,LR.X1.X2=lm2))
summary(r1)
Brier(list(lm1,lm2))

# machine learning
library(randomForest)
dat$Y=factor(dat$Y)
rf <- randomForest(Y~X2,data=dat)
rocCV=Roc(list(RandomForest=rf,LogisticRegression=lm2),
    data=dat,
    splitMethod="bootcv",
    B=3,
    cbRatio=1)
plot(rocCV)

# compute .632+ estimate of Brier score
bs <- Brier(list(LR.X1=lm1,LR.X2=lm2),
    data=dat,
    splitMethod="boot632+",
    B=3)
bs
#'
}
\author{
Thomas Gerds \email{tag@biostat.ku.dk}
}
\references{
Fawcett, T. (2006). An introduction to ROC analysis.
Pattern Recognition Letters, 27, 861-874.

Gerds, Cai & Schumacher (2008). The Performance of Risk
Prediction Models. Biometrical Journal, Vol 50, 4, 457-479.

Efron, Tibshirani (1997) Journal of the American
Statistical Association 92, 548--560 Improvement On
Cross-Validation: The .632+ Bootstrap Method.

Wehberg, S and Schumacher, M (2004) A comparison of
nonparametric error rate estimation methods in
classification problems. Biometrical Journal, Vol 46,
35--47
}
\keyword{models}

