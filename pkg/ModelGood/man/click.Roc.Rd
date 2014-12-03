\name{click.Roc}
\alias{click.Roc}
\title{Click on ROC curve}
\usage{
click.Roc(object, pch = 19, label = TRUE, adj, col = "orange", cex = 3,
  ...)
}
\arguments{
  \item{object}{An object obtained with function
  \code{Roc}}

  \item{pch}{the symbol of the tag}

  \item{label}{If TRUE label the tag.}

  \item{adj}{passed to \code{text} to adjust of the legend
  relative to clickpoint.}

  \item{col}{the color of the tag}

  \item{cex}{the size of the tag}

  \item{...}{passed to \code{identify}}
}
\value{
the values at the tag
}
\description{
Show marker value sensitivity and specificity at mouse
point
}
\details{
A tag is set on the ROC curve at the mouse click and
corresponding marker value, sensitivity and specificity
shown below the click-point.
}
\examples{
\dontrun{
x <- abs(rnorm(20))
d <- data.frame(y=rbinom(1:20,1,p=x/max(x)))
r <- Roc(y~x,data=d)
plot(r)
click.Roc(r)
 }
}
\author{
Thomas A. Gerds <tag@biostat.ku.dk>
}
\seealso{
identify Roc
}

