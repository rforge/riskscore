MgResolveSmartArgs <- function(call,
                               keys,
                               ignore,
                               defaults,
                               forced,
                               split="\\.",
                               verbose=TRUE){
  SmartArgs <- as.list(call)
  # remove ignorable arguments
  SmartArgs <- SmartArgs[names(SmartArgs)!=""]
  # ignore case
  names(SmartArgs) <- tolower(names(SmartArgs))
  if (!missing(ignore) && is.character(ignore))
    SmartArgs <- SmartArgs[match(names(SmartArgs),ignore,nomatch=0)==0]
  if (verbose==TRUE){
    allKeysRegexp <- paste("^",keys,split,sep="",collapse="|")
    notIgnored <- grep(allKeysRegexp,names(SmartArgs),value=FALSE,ignore.case=TRUE)
    Ignored <- names(SmartArgs)[-notIgnored]
    SmartArgs <- SmartArgs[notIgnored]
    if (length(Ignored)>0)
      warning(paste("The following argument(s) are not smart and therefore ignored: ",paste(Ignored,collapse=", ")))
  }
  # ---------------------------default arguments------------------------
  keyDefaults <- vector(mode="list",length=length(keys))
  names(keyDefaults) <- keys
  if (!missing(defaults)){
    whereDefault <- match(names(defaults),names(keyDefaults),nomatch=0)
    if (all(whereDefault))
      keyDefaults[whereDefault] <- defaults
    else stop("Not all default arguments found.")
  }
  # ---------------------------forced arguments------------------------
  keyForced <- vector(mode="list",length=length(keys))
  names(keyForced) <- keys
  if (!missing(forced)){
    whereDefault <- match(names(forced),names(keyForced),nomatch=0)
    if (all(whereDefault))
      keyForced[whereDefault] <- forced
    else stop("Not all forced arguments found.")
  }
  # ------------------extract args that match key------------------
  keyArgList <- lapply(keys,function(k){
    keyRegexp <- paste("^",k,split,sep="")
    foundArgs <- grep(keyRegexp,names(SmartArgs),value=TRUE,ignore.case=TRUE)
    if (length(foundArgs)>0){
      keyArgs <- SmartArgs[foundArgs]
      argNames <- sapply(strsplit(names(keyArgs),keyRegexp),function(x)x[[2]])
      keyArgs <- lapply(keyArgs,eval)
      names(keyArgs) <- argNames
    }
    else{
      keyArgs <- NULL
    }
    # -----------------prepending the forced arguments-----------------

    if (length(keyForced[[k]])>0){
      keyArgs <- c(keyForced[[k]],keyArgs)
    }
    # -----------------appending the default arguments-----------------
    if (length(keyDefaults[[k]])>0){
      keyArgs <- c(keyArgs,keyDefaults[[k]])
    }
    # ------------------------removing duplicates------------------------
    keyArgs[!duplicated(names(keyArgs))]
  })

  names(keyArgList) <- keys
  keyArgList
}

## dummy <- function(x,...){
##   resolveSmartArgs(call=match.call(),keys=c("test","surv","aa.bb"),defaults=list("test"=list(u=7)),forced=list("test"=list(v=1),"surv"=list(x=1)),ignore="x")
## }
## dummy(x=1,test.v=2,aa.bb.cc=22)
## dummy(x=1,test.V=2,aa.bb.cc=22)
## dummy(x=1,aa.bb.cc=22)
