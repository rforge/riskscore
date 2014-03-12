"summary.Brier" <- function(object,digits=2,print.it=TRUE,...){
  if (!is.null(object$method) && print.it) print(object$method)
  res <- c("BS","AppBS","BootcvBS","NoInfBS")
  names(res) <- c(object$method$name,c("App","BCV","NoInf"))
  found <- match(names(object),res,nomatch=FALSE)
  res <- res[found]
  out <- lapply(res,function(r){
    out <- sapply(1:length(object$models),function(w){
      object[[r]][[w]]
    })
    names(out) <- names(object$models)
    out
  })
  outMat <- 100*do.call("cbind",out)
  outMat <- outMat[order(outMat[,NCOL(outMat)]),,drop=FALSE]
  if (print.it) cat("\nEstimated Brier score in %\n")
  if (print.it) print(apply(as.data.frame(outMat),2,round,digits=digits),quote=FALSE)
  invisible(out)
}
