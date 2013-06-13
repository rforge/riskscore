library(COST)
library(pec) 
library(randomForest)
library(randomForestSRC)
library(riskRegression) 
library(cmprsk)
source("~/Dropbox/PseudoForests/R/pseudoForest.R")
source("~/Dropbox/PseudoForests/R/predictEventProb.pseudoForest.R")
data(cost,package="COST")
cost$cause=ifelse(cost$event %in% c(1,2,3,6), 1, ifelse(cost$event!=0,2,0)) 
cost0 <- na.omit(cost)
fixHorizon=5*365.25

## source("~/research/riskPlots/R/predictStatusProb.R")
## source("~/research/riskPlots/R/riskplot.R")
## source("~/research/riskPlots/R/riskplot.R")

# {{{ Uncensored

logreg = glm(formula=am ~ hp + wt,data=mtcars,family=binomial) 

## predictStatusProb(logreg,newdata=mtcars[1:2,])
## newD <- data.frame(wt=1.7,hp=200)
## predictStatusProb(logreg,newdata=newD)
## predict(logreg,newdata=mtcars[1:2,],type="response")

## source("~/research/riskPlots/R/riskplot.R")
riskplot(logreg,am ~ hp + wt,data=mtcars,plotObs=TRUE)

## Old version riskplot
## source("~/research/riskPlots/R/riskplot.R")
## riskplot(logreg,am ~ inX(hp) + inY(wt),data=mtcars)
## riskplot(am ~ inX(hp)+inY(wt),logreg)


# }}}
# {{{ Survival example

# survival example
survForm = Surv(time,status) ~ age + strokeScore
cox <- cph(survForm, data=cost0,surv=TRUE)
rsfSurv <- rfsrc(Surv(time,status) ~ age + strokeScore, data=cost0,ntree=50,keep.forest=TRUE)

#COX
## source("~/research/riskPlots/R/riskplot.R")
riskplot(cox,Hist(time,status)~strokeScore+age,data=cost0,horizon=fixHorizon)

## FIX
## riskplot(cox,Hist(time,status)~time+age,data=cost0,horizon=fixHorizon)

#RSF
riskplot(rsfSurv,Hist(time,status)~strokeScore+age,data=cost0,horizon=fixHorizon)


# }}}
# {{{ Competing risks example

# Competing risk example
cpForm=Hist(time,cause)~age + strokeScore
fixCause=1
fitprf <- pseudoForest(cpForm,data=cost0,cause=fixCause,times=c(0,fixHorizon),ntree=1000,keep.forest=TRUE,replace=FALSE)
fitcsc <- CSC(cpForm, data=cost0, survtype="hazard")

riskplot(fitcsc,
             Hist(time,cause)~strokeScore+age,
             data=cost0,
             horizon=fixHorizon,
             cause=1)


riskplot(fitprf,
             Hist(time,cause)~strokeScore+age,
             data=cost0,
             horizon=fixHorizon,
             cause=2)


# }}}

# {{{ linear predictor

aa <- attr(terms(formula(logreg$call)),"term.labels")
da <- eval(logreg$call$data)[,aa]

predictStatusProb(logreg,newdata=da)
rowSums(predict(logreg,type="terms"))

linPredObs <- function(object){
  beta <- coefficients(object)
  
  beta
}

# }}}


# --------------------------ORIGINAL version--------------------------
## source("~/research/riskPlots/R/riskplot.R")
## riskplotForm=~inY(age)+inX(strokeScore)
## riskplot(object=fitprf,formula=riskplotForm,ref="1",control=list(colorkey=TRUE),cause=1,xlab="Scandinavian stroke score",ylab="Age (years)",main="Pseudo random forest",horizon=fixHorizon,response=with(cost0,Hist(time,cause)),pchObs=as.numeric(as.character(factor(with(cost0,Hist(time,cause))[,"event"],levels=1:3,labels=c(16,1,2),))),data=cost0) 
