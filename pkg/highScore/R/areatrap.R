areatrap<-function(Abs,Ord){
  nobs <- length(Abs)
  dAbs <- Abs[-1]-Abs[-nobs]
  mil <- (Ord[-nobs]+Ord[-1])/2
  area <- sum(dAbs*mil)
  return(area)
}
