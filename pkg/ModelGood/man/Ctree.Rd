\name{Ctree}
\alias{Ctree}
\title{S3 wrapper function for party's ctree method}
\usage{
Ctree(...)
}
\arguments{
  \item{...}{passed to \code{ctree}}
}
\value{
object of class Ctree which contains a ctree object
}
\description{
S3 wrapper function for party's ctree method
}
\details{
The ModelGood crossvalidation functionality works only for
S3 classes.
}
\examples{
library(party)
set.seed(7)
x <- abs(rnorm(20))
d <- data.frame(y=rbinom(20,1,x/max(x)),x=x,z=rnorm(20))
ct <- Ctree(y~x+z,d)
plot(ct$ctree)
predictStatusProb(ct,newdata=d[1:3,])
}
\author{
Thomas A. Gerds <tag@biostat.ku.dk>
}
\seealso{
ctree
}

