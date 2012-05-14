# {{{ Response

library(survival)
library(pec)
check.code("highScore")
set.seed(177)
d <- SimSurv(10)
d2 <- prodlim:::SimCompRisk(10)

f <- coxph(Surv(time,status)~X1,data=d)

score(object=list(f),formula=time~1,data=d)
score(object=list(f),formula=status~1,data=d)
d2$dummy <- factor(d2$cause)
score(object=list(f),formula=dummy~1,data=d2)

score(object=list(f),formula=Surv(time,status)~1,data=d)
score(object=list(f),formula=Surv(time,status)~1,data=d,censMethod="pseudo")
score(object=list(f),formula=Surv(time,status)~1,data=d[d$status!=0,],censMethod="pseudo")

score(object=list(f),formula=Hist(time,cause)~1,data=d2)

h2 <- with(d2,Hist(time,cause))

score(object=list(f))

# }}}


