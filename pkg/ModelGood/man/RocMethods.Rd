\name{RocBrierScore}
\alias{Brier}
\alias{Brier.list}
\alias{Brier.glm}
\alias{Brier.lrm}
\alias{Brier.randomForest}
\alias{Brier.rpart}
\alias{Roc}
\alias{Roc.list}
\alias{Roc.glm}
\alias{Roc.lrm}
\alias{Roc.randomForest}
\alias{Roc.rpart}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Comparing prediction models with Receiver operating characteristics and Brier's scores}
\description{
  Evaluating the performance of risk prediction models with binary
  status response variable (case/control or similar). Roc curves are based on a single marker,
  which is the probability predicted case probability by a model. The 
   area under the curve and the Brier score is used to summarize the performance.

  Bootstrap-crossvalidation techniques are implemented to
  estimate the generalization performance of the model(s), i.e. the performance which
  can be expected in new subjects.
}

\usage{
Roc(object,...)
\method{Roc}{list}(object,
                     formula,
                     data,
                     splitMethod="noSplitMethod",
                     noinf.method=c("simulate"),
                     simulate="reeval",
                     B,
                     M,
                     breaks,
                     crRatio=1,
                     RocAverageMethod="vertical",
                     RocAverageGrid=switch(RocAverageMethod,
		       "vertical"=seq(0,1,.01),
		       "horizontal"=seq(1,0,-.01)),
                     model.args=NULL,
                     model.parms=NULL,
                     keepModels=FALSE,
                     keepSampleIndex=FALSE,
                     keepCrossValRes=FALSE,
                     keepNoInfSimu,
                     slaveseed,
                     na.accept=0,
                     verbose=TRUE,
                     ...)
Brier(object,...)
\method{Brier}{list}(object,
                     formula,
                     data,
                     splitMethod="noSplitMethod",
                     noinf.method=c("simulate"),
                     simulate="reeval",
                     crRatio=1,
                     B,
                     M,
                     model.args=NULL,
                     model.parms=NULL,
                     keepModels=FALSE,
                     keepSampleIndex=FALSE,
                     keepCrossValRes=FALSE,
                     na.accept=0,
                     verbose=TRUE,
                     ...)
\method{Brier}{glm}(object,formula,data,...)
\method{Brier}{lrm}(object,formula,data,...)
\method{Brier}{rpart}(object,formula,data,...)
\method{Brier}{randomForest}(object,formula,data,...)
\method{Roc}{glm}(object,formula,data,...)
\method{Roc}{lrm}(object,formula,data,...)
\method{Roc}{rpart}(object,formula,data,...)
\method{Roc}{randomForest}(object,formula,data,...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A named list of prediction models. Each entry
    is either an R-object for which a \link{predictStatusProb}
    method exists (see details) or a \code{call} that can be evaluated to such an R-object.
    For cross-validation (e.g. when \code{splitMethod} is "boot632plus")
    all the models in this list must include their \code{call} in their
    value.
  }
  \item{formula}{
    A formula. If missing, use the formula of the
    (first) model in object.
    The left hand side is used to find the status response
    variable in \code{data}.
  }
  \item{data}{A data frame to validate the prediction models
    If missing, use the data of the (first) model in object.
  }
  \item{splitMethod}{Method for estimating the generalization error.

    \code{none}:Assess the models in the same data where they are
      fitted. Yields the apparent or re-substitution performance. Often
      overestimates the generalization performance.

      \code{BootCV}: Bootstrap cross validation. The prediction models are
      trained on \code{B} bootstrap samples, that are either drawn with or
      without replacement from \code{data} of the size \code{M}.
      The model performance is estimated with  the observations that are NOT in the
      bootstrap sample.

    \code{Boot632}: Linear combination of the apparent performance and
      the  Bootcv performance using
      the constant weight .632 (see references).

    \code{Boot632plus}: Linear combination of apparent performance and Bootcv using
    weights dependent on how the models perform in permuted data (see
    references). 
      
    \code{NoInf}: 
    Assess the models in permuted data.
    }
  \item{noinf.method}{Method for computing the no-information Roc
    curve.}
  \item{simulate}{If equal to \code{"reeval"} then the models are re-evaluated in the current
   permuted data for computing the no-information Roc curve.}
  \item{B}{Number of bootstrap samples. The default depends on argument \code{splitMethod}.
    When \code{splitMethod in c("Bootcv","Boot632","Boot632plus"} the
    default is 100. For \code{splitMethod="cvK"} \code{B} is the number of cross-validation cycles,
    and -- default is 1.
    For \code{splitMethod="none"} \code{B} is the number of bootstrap
    simulations e.g. to obtain bootstrap confidence limits -- default is 0.}
  \item{M}{The size of the bootstrap samples for cross-validation without replacement.
  }
  \item{breaks}{Break points for computing the Roc curve. Defaults to
    \code{seq(0,1,.01)} for the Roc.list method and to
    \code{sort(unique(breaks))} for the default method.
  }
  \item{crRatio}{Experimental and not yet tested. Cost/benefit ratio. Default value is to 1, meaning that
    misclassified cases are as bad as misclassified controls.}
  \item{RocAverageMethod}{Method for averaging the bootstrapped Roc
    curves. Either \code{"vertical"} or
    \code{"horizontal"}. See reference 1 below for details.}
  \item{RocAverageGrid}{Grid points for the averaging of boostrapped Roc
  curves.}
  \item{model.args}{List of extra arguments that
    can be passed to the \code{predictStatusProb} methods. The list must have
    an entry for each entry in \code{object}.
  }
  \item{model.parms}{
    List with exactly one entry for each entry in \code{object}.
    Each entry names parts of the value of the fitted models that should
    be extracted and added to the output (see value). 
  }
  \item{keepModels}{If \code{FALSE} keep only the names of the elements
    of object.
    If \code{"Call"} then keep the call of the elements of object.
    Else, add the object as it is to the output.
  }
  \item{keepSampleIndex}{
    Logical. If \code{FALSE} remove the cross-validation
    index (which tells who was in the learn and who in the validation
    set)
    from the output list which otherwise is included in the method
    part of the output list. 
  }
  \item{keepCrossValRes}{
    Logical. If \code{TRUE} add all \code{B} crossvalidation results to
     the output (see value). Defaults to \code{TRUE}.
  }
  \item{keepNoInfSimu}{
    Logical. If \code{TRUE} add the \code{B} results in permuted data
    (for no-information performance) to
     the output (see value). Defaults to \code{FALSE}.
  }
  \item{slaveseed}{
    Vector of seeds, as long as \code{B}, to be given to the slaves in parallel computing.
  }
  \item{na.accept}{
    Works only for "Bootcv" estimate of performance.
    The maximal number of failures during training
    the models to the bootstrap samples. Usually a
    small number relative to \code{B}.
  }
  \item{verbose}{if \code{TRUE} the procedure is reporting details of
    the progress, e.g. it prints the current step in cross-validation
    procedures.
  }
  \item{...}{Difficult to explain}
}
\details{
  Missing data in the response or in the input matrix cause a failure. 
  
  For each prediction model there must
  be a \code{predictStatusProb} method: for example, to assess a 
  prediction model which evaluates to a \code{myclass} object one defines 
  a function called \code{predictStatusProb.myclass} with arguments
  \code{object,newdata,cutpoints,train.data,...}, like this

  myFit=myModel(Y~X,data=dat)
  
  class(myFit)="myclass"

  \code{predictStatusProb.myclass <-
  function(object,newdata,cutpoints,train.data,...){
    predict(object, data=newdata,method="probabilities")
    out
  }}

  Such a function takes the object which was fitted with train.data and
  derives a matrix with predicted event status probabilities for each subject
  in newdata (rows) and each cutpoint (column) of the response variable
  that defines an event status.
  
  Currently implemented are \code{predictStatusProb} methods  for the following R-functions:
  \describe{
    \item{}{\code{glm} (from \code{library(stats)}}
    \item{}{\code{lrm} (from \code{library(Design)}}
    \item{}{\code{rpart} (from \code{library(rpart)})}
    \item{}{\code{randomForest} from \code{library(randomForest)}}
  }
}

\value{
  Object of class \code{Roc} or class \code{Brier} for which 
   \code{print.Roc}, \code{summary.Roc},
  and \code{plot.Roc} (only \code{Roc}) methods are available.
  
  The object includes the following components:
  \item{Roc}{
    The Roc curve(s) estimated according to the
    \code{splitMethod}. A matrix where each column represents the
    estimated prediction error of a fit at
    the time points in time.
  }
  \item{AppRoc}{
    The Roc curve(s) estimated when
    the model(s) are evaluated in the same data where
    they were fitted. Only if \code{splitMethod} is one of
    "NoInf", "Bootcv", "Boot632" or "Boot632plus", since otherwise
    \code{repla} is "apparent" and then this is 
     stored in \code{Roc} as explained above.
  }
  \item{BootcvRoc}{
    The prediction error when
    the model(s) are trained in the bootstrap sample and evaluated in
    the data that are not in the bootstrap sample.
    Only if \code{splitMethod} is one of
    "Boot632" or "Boot632plus". When \code{splitMethod="Bootcv"}
    then the \code{BootcvRoc} is stored in the component \code{PredRoc}.
  }
  \item{NoInfRoc}{
    The prediction error when
    the model(s) are evaluated in the permuted data.
    Only if \code{splitMethod} is one of
    "Bootcv", "Boot632", or "Boot632plus".
    For \code{splitMethod="NoInf"}
    the \code{NoInfRoc} is stored in the component \code{PredRoc}.
  }
  \item{weight}{
    The weight used to linear combine the \code{AppRoc} and the
    \code{BootcvRoc} 
    Only if \code{splitMethod} is one of
    "Boot632", or "Boot632plus".
  }
  \item{overfit}{
    Estimated \code{overfit} of the model(s).
    See references.
    Only if \code{splitMethod} is one of
    "Boot632", or "Boot632plus".
  }
  \item{call}{The call that produced the object} 
  \item{models}{See keepModels}
  \item{method}{The method used for estimation of the overfitting bias.}
}

\references{
  Fawcett, T. (2006). An introduction to ROC analysis. Pattern
  Recognition Letters, 27, 861-874.
  
  Gerds, Cai & Schumacher (2008). The Performance of Risk Prediction
  Models. Biometrical Journal, Vol 50, 4, 457-479.

  Efron, Tibshirani (1997) Journal of the American Statistical Association 92, 548--560
  Improvement On Cross-Validation: The .632+ Bootstrap Method.

  
  Wehberg, S and Schumacher, M (2004)
  A comparison of nonparametric error rate estimation methods in classification problems
  Biometrical Journal, Vol 46,  35--47
  }

\author{Thomas Gerds \email{tag@biostat.ku.dk} }

% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{} % plot.Roc, \code{\link{summary.Roc}}}
\examples{
## Generate some data with binary response Y
## depending on X1 and X2 and X1*X2 
N <- 400
X1 <- rnorm(N)
X2 <- rbinom(N,1,.4)
expit <- function(x) exp(x)/(1+exp(x))
lp <- expit(1 + X1 + X2 - X1*X2)
Y <- factor(rbinom(N,1,lp))
dat <- data.frame(Y=Y,X1=X1,X2=X2)

## fit a logistic model
lm1 <- glm(Y~X1,data=dat,family="binomial")
lm2 <- glm(Y~X1+X2,data=dat,family="binomial")
r1=Roc(list(lm1,lm2),verbose=0,crRatio=1)
summary(r1)
Brier(list(lm1,lm2),verbose=0,crRatio=1)


# crossing curves
set.seed(18)
N=500
Y=rbinom(N,1,.5)
X1=rnorm(N)
X1[Y==1]=rnorm(sum(Y==1),mean=rbinom(sum(Y==1),1,.5))
X2=rnorm(N)
X2[Y==0]=rnorm(sum(Y==0),mean=rbinom(sum(Y==0),1,.5))
dat <- data.frame(Y=Y,X1=X1,X2=X2)
lm1 <- glm(Y~X1,data=dat,family="binomial")
lm2 <- glm(Y~X2,data=dat,family="binomial")
plot(Roc(list(lm1,lm2),data=dat,verbose=0,crRatio=1))
\dontrun{
library(randomForest)
dat$Y=factor(dat$Y)
rf <- randomForest(Y~X2,data=dat)
rocCV=Roc(list(RandomForest=rf,LogisticRegression=lm2),
  data=dat,
  verbose=TRUE,
  splitMethod="bootcv",
  B=10,
  crRatio=1)
plot(rocCV,diag=TRUE)
}
 }
\keyword{models}% at least one, from doc/KEYWORDS
