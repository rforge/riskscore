\name{selectFGR}
\alias{selectFGR}
\title{
Stepwise variable selection in the Fine & Gray regression competing risk model 
}
\description{
  This is a wrapper function which first selects variables in the Fine & Gray 
  regression model using \code{crrstep} from the \code{crrstep} package and 
  then returns a fitted Fine & Gray regression model with the selected variables.
}
\usage{
selectFGR(formula, data, cause=1, rule = "AIC", 
          direction="backward", ...)
}
\arguments{
    \item{formula}{
    A formula whose left hand side is 
    a \code{Hist} object -- see \code{\link{Hist}}.
    The right hand side specifies (a linear combination of)
    the covariates. See examples below.
  }
    \item{data}{
    A data.frame in which all the variables of
    \code{formula} can be interpreted.
  }
  \item{cause}{
    The failure type of interest. Defaults to \code{1}. 
  }
  \item{rule}{Rule to pass on to crrstep ("AIC", "BIC" or "BICcr"), also
    see \code{crrstep}}
  \item{direction}{see \code{crrstep}}
  \item{\dots}{
    Further arguments passed to \code{crrstep}. 
}
}
\author{
  Rob C.M. van Kruijsdijk \email{R.C.M.vanKruijsdijk@umcutrecht.nl}
  
  Thomas Alexander Gerds \email{tag@biostat.ku.dk}
}
\examples{
\dontrun{
library(riskRegression)
library(prodlim)
library(pec)
library(cmprsk)

d <- SimCompRisk(100)
newvars <- as.data.frame(matrix(runif(8*100),nrow=100))
colnames(newvars) <- c('X3','X4','X5','X6','X7','X8','X9','X10')
d <- cbind(d,newvars)
fg <- FGR(Hist(time, cause) ~ X1 + X2 + X3 + X4,cause=1,data=d)
sfg <- selectFGR(Hist(time, cause) ~ X1 + X2 + X3 + X4,cause=1,data=d,rule="AIC",direction="backward") 
predictEventProb(sfg,times=100,newdata=d)

## apparent performance
pec(list(full.model=fg,selected.model=sfg),formula=Hist(time, cause)~1,data=d)

## bootstrap cross-validation performance
set.seed(7)
pec(list(full.model=fg,selected.model=sfg),formula=Hist(time, cause)~1,data=d,B=10,splitMethod="bootcv")
}
}

\keyword{survival}
