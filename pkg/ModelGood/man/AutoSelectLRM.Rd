\name{AutoSelectLRM}
\alias{AutoSelectLRM}
\title{Automated backward elemination for logistic regression}
\usage{
AutoSelectLRM(formula, data, ...)
}
\arguments{
  \item{formula}{passed to lrm}

  \item{data}{passed to lrm}

  \item{...}{passed to fastbw}
}
\value{
object of class AutoSelectLRM
}
\description{
Wrapper for automated backward elemination for logistic
regression
}
\details{
First run backward elemination via fastbw from the rms
package, then fit the logistic regression model including
the selected variables
}
\examples{
library(rms)
set.seed(7)
x <- abs(rnorm(20))
d <- data.frame(y=rbinom(20,1,x/max(x)),x=x,z=rnorm(20))
fbw <- AutoSelectLRM(y~x+z,d)
predictStatusProb(fbw,newdata=d[1:3,])
}
\author{
Thomas A. Gerds <tag@biostat.ku.dk>
}
\seealso{
fastbw lrm
}

