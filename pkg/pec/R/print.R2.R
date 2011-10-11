
print.R2 <- function(x,...){
  cat("\nTime-dependent explained variation:\n\n 1- Brier(model)/Brier(reference)\n\nReference model: ",attr(x,"reference"),"\n\n")
  ## if (all(sapply(x,length)==1)){
  ## x <- do.call("rbind",x)
  ## }
  print.listof(x,...)
}
