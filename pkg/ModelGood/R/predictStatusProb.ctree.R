# S3 wrapper function for party's ctree method
# --------------------------------------------------------------------
MGctree <- function(...){
    ## require(party)
    ctree <- party::ctree
    out <- list(ctree=ctree(...))
    class(out) <- "MGctree"
    out$call <- match.call()
    out  
}

predictStatusProb.MGctree <- function (object, newdata, ...) {
 ## require(party)
 N <- NROW(newdata)
 p <- party::treeresponse(object$ctree, newdata=newdata)
 if (NROW(p) != NROW(newdata))
   stop("Prediction failed")
 p
}
