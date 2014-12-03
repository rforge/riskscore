\name{ElasticNet}
\alias{ElasticNet}
\title{Wrapper function for glmnet}
\usage{
ElasticNet(formula, data, nfolds = 10, ...)
}
\arguments{
  \item{formula}{Formula where the right hand side
  specifies the response and the left hand side the
  predictor matrix}

  \item{data}{A data frame in which \code{formula} is
  evaluated}

  \item{nfolds}{nfolds: number of cross-validation folds in
  cv.glmnet (default in function is 10)}

  \item{...}{passed on to glmnet}
}
\value{
Object with class ElasticNet
}
\description{
Wrapper function for glmnet
}
\details{
This function first calls \code{cv.glmnet} and then
evaluates glmnet at the hyper parameter which optimizes the
cross-validation criterion.
}
\examples{
# Generate some data with binary response Y
 # depending on X1 and X2 and X1*X2
set.seed(40)
N <- 40
X1 <- rnorm(N)
X2 <- rbinom(N,1,.4)
X3 <- rnorm(N)
expit <- function(x) exp(x)/(1+exp(x))
lp <- expit(1 + X1 + X2 + X3)
Y <- factor(rbinom(N,1,lp))
dat <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)

efit <- ElasticNet(Y~X1+X2+X3,data=dat,family="binomial",alpha=0.1)
Brier(efit,verbose=FALSE)
}
\author{
Thomas A. Gerds <tag@biostat.ku.dk>
}
\seealso{
predictStatusProb
}

