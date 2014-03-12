"MgSindex2" <- function (jump.times, eval.times, comp="smaller",strict=FALSE) {
  stopifnot(is.numeric(jump.times))
  stopifnot(is.numeric(eval.times))
  N <- length(jump.times)
  if (comp=="greater"){
    N-MgSindex(jump.times=jump.times,eval.times=eval.times,comp="smaller",strict=!strict)
  }
  else{
    neval <- length(eval.times)
    if (!(neval> 0 && N >0)) stop("missing data")
    new.order <- order(eval.times)
    ind <- .C("MgSindex",index = integer(neval),as.double(sort(jump.times)),as.double(eval.times[new.order]),as.integer(N),as.integer(neval),as.integer(strict),DUP=FALSE,PACKAGE="ModelGood")$index
    ind[order(new.order)]
  }
}
