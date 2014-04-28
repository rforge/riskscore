\docType{data}
\name{Pbc3}
\alias{Pbc3}
\title{Pbc3 data}
\format{A data frame with 349 observations on the following 15 variables.
\describe{ \item{list("ptno")}{patient identification}
\item{list("unit")}{hospital (1: Hvidovre, 2: London, 3: Copenhagen, 4:
Barcelona, 5: Munich, 6: Lyon)} \item{list("tment")}{treatment (0: placebo,
1: CyA)} \item{list("sex")}{(1: males, 0: females)} \item{list("age")}{age
in years} \item{list("stage")}{histological stage (1, 2, 3, 4)}
\item{list("gibleed")}{previous gastrointestinal bleeding (1: yes, 0: no)}
\item{list("crea")}{creatinine (micromoles/L)} \item{list("alb")}{albumin
(g/L)} \item{list("bili")}{bilirubin (micromoles/L)}
\item{list("alkph")}{alkaline phosphatase (IU/L)}
\item{list("asptr")}{aspartate transaminase (IU/L)}
\item{list("weight")}{body weight (kg)} \item{list("days")}{observation time
(days)} \item{list("status")}{status at observation time (0: censored, 1:
liver transplantation, 2 : dead)} }}
\source{
\url{http://192.38.117.59/~linearpredictors/?page=datasets&dataset=Pbc3}
}
\description{
PBC3 was a multi-centre randomized clinical trial conducted
in six European hospitals. Between 1 Jan. 1983 and 1 Jan.
1987, 349 patients with the liver disease primary biliary
cirrhosis (PBC) were randomized to either treatment with
Cyclosporin A (CyA, 176 patients) or placebo (173
patients). The purpose of the trial was to study the effect
of treatment on the survival time. However, during the
course of the trial an increased use of liver
transplantation for patients with this disease made the
investigators redefine the main response variable to be
time to ``failure of medical treatment'' defined as either
death or liver transplantation. Patients were then followed
from randomization until treatment failure, drop-out or 1
Jan, 1989; 61 patients died (CyA: 30, placebo: 31), another
29 were transplanted (CyA: 14, placebo: 15) and 4 patients
were lost to follow-up before 1 Jan. 1989. At entry a
number of clinical, biochemical and histological variables,
including serum bilirubin, serum albumin, sex, age were
recorded.
}
\examples{
data(Pbc3)
}
\references{
Andersen and Skovgaard. Regression with linear predictors.
Springer, 2010.
}
\keyword{datasets}

