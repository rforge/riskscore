FastBw <- function(formula,data,rule="aic"){
  fit <- lrm(formula,data)
  bwfit <- fastbw(fit,rule=rule)
  if (length(bwfit$names.kept)==0){
    newform <- reformulate("1",formula[[2]])
    newfit <- glm(newform,data,family="binomial")
  }
  else{
    newform <- reformulate(bwfit$names.kept, formula[[2]])
    newfit <- lrm(newform,data)
  }
  out <- list(fit=newfit,In=bwfit$names.kept)
  out$call <- match.call()
  class(out) <- "FastBw"
  out
}
predictStatusProb.FastBw <- function(object,newdata,...){
  predictStatusProb(object[[1]],newdata=newdata,...)
}
