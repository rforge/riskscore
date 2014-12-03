\name{Sensitivity}
\alias{Diagnose}
\alias{NPV}
\alias{PPV}
\alias{Sensitivity}
\alias{Specificity}
\title{Compute sensitivity, specificity and predictive values}
\usage{
Sensitivity(x,event,cutoff,comparison=">=",...)
       Specificity(x,event,cutoff,comparison=">=",...)
       NPV(x,event,cutoff,comparison=">=",...)
       PPV(x,event,cutoff,comparison=">=",...)
}
\arguments{
  \item{x}{Either a binary 0,1 variable, or a numeric
  marker which is cut into binary.}

  \item{event}{Binary response variable. Either a 0,1
  variable where 1 means 'event', or a factor where the
  second level means 'event'.}

  \item{cutoff}{When x is a numeric marker, it is compared
  to this cutoff to obtain a binary test.}

  \item{comparison}{How x is to be compared to the cutoff
  value}

  \item{...}{passed on to \code{binom.test}}
}
\value{
list with Sensitivity, Specificity, NPV, PPV and confidence
interval
}
\description{
Compute sensitivity, specificity and predictive values
}
\details{
Confidence intervals are obtained with \code{binom.test}
}
\examples{
set.seed(17)
x <- rnorm(10)
y <- rbinom(10,1,0.4)
Sensitivity(x,y,0.3)
Specificity(x,y,0.3)
PPV(x,y,0.3)
NPV(x,y,0.3)

Diagnose(x,y,0.3)
}
\author{
Thomas A. Gerds <tag@biostat.ku.dk>
}
\seealso{
binom.test
}

