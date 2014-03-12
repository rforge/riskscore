summary.Roc <- function(object,digits=2,print.it=TRUE,...){
  if (print.it) cat("\nReceiver operating characteristic\n")
  if (!is.null(object$method) && print.it) print(object$method)
  if (print.it){cat("\nResponse:\n")
                print(table(object$Response))
                cat("\n")}
  if (length(object$models)>0){
    
    # ------------------------Area under the curve------------------------

    modelsAuc <- do.call("rbind",lapply(object$Auc,function(m){round(100*unlist(m),digits=digits)}))
    rownames(modelsAuc) <- names(object$models)
    xxx <- names(object$Auc[[1]])
    names(xxx) <- xxx
    xxx[xxx=="Auc"] <- paste("AUC",object$method$name,sep=":")
    colnames(modelsAuc) <- xxx
    ##     colnames(modelsAuc) <- names(object$Auc[[1]])
    if (print.it) cat("\nArea under the curve:\n\n")
    ## print(apply(as.data.frame(modelsAuc[order(-modelsAuc[,1,drop=TRUE]),,drop=FALSE]),2,formatC,digits=digits),quote=FALSE)
    if (print.it)
      print(apply(as.data.frame(modelsAuc[order(-modelsAuc[,1,drop=TRUE]),,drop=FALSE]),2,round,digits=digits),quote=FALSE)
    
    
    # ----------------------------Brier score----------------------------

    modelsBS <- do.call("rbind",lapply(object$Brier,function(m){round(100*unlist(m),digits=digits)}))
    rownames(modelsBS) <- names(object$models)
    ## colnames(modelsBS) <- c(paste("BS",object$method$name,sep=":"),names(object$Brier[[1]])[-1])
    ## colnames(modelsBS) <- names(object$Brier[[1]])
    xxx <- names(object$Brier[[1]])
    names(xxx) <- xxx
    xxx[xxx=="BS"] <- paste("BS",object$method$name,sep=":")
    colnames(modelsBS) <- xxx
    if (print.it) cat("\n\nBrier score:\n\n")
    ## print(apply(as.data.frame(modelsBS[order(modelsBS[,1,drop=TRUE]),,drop=FALSE]),2,formatC,digits=digits),quote=FALSE)
    if (print.it)
      print(apply(as.data.frame(modelsBS[order(modelsBS[,1,drop=TRUE]),,drop=FALSE]),2,round,digits=digits),
            quote=FALSE)
    out <- list(Auc=modelsAuc,Brier=modelsBS)
    class(out) <- "summary.Roc"
    invisible(out)
  }
  else{
    stop("Dont know what happened.")
    ## cat(round(100*Auc.default(object$Roc$Sensitivity,object$Roc$Specificity),digits=digits))
  }
}

print.summary.Roc <- function(x,digits=2,...){
  lapply(x,function(x){
    cat("\n")
    print(apply(as.data.frame(x),2,round,digits=digits),quote=FALSE)
  })
}

##   cat("\n\n")
##   if (!is.null(object$PPV)){
##     PV(x)


## PV <- function(object)
## {
##   if (!is.null(object$confint) && object$confint==TRUE){
##     cat("Confidence method: ",object$confint.method)
##     cat("\n\n")
##     prSens <- object$Roc$Sensitivity[who]
##     prSpec <- object$Roc$Specificity[who]
##     prPPV <- object$Roc$PPV[who]
##     prNPV <- object$Roc$NPV[who]
##     myDigits <- if (percent==TRUE) digits-2 else digits
##     myFactor <- if (percent==TRUE) 100 else 1 
##     out <- do.call("cbind",list(Sens=round(myFactor*prSens,myDigits),"(CI.Sens)"=MgFormCi(lower=object$CI.Sens[who,"Lower"],upper=object$CI.Sens[who,"Upper"]),Spec=round(myFactor*prSpec,myDigits),"(CI.Spec)"=MgFormCi(lower=object$CI.Spec[who,"Lower"],upper=object$CI.Spec[who,"Upper"]),PPV=round(myFactor*prPPV,myDigits),"(CI.PPV)"=MgFormCi(lower=object$CI.PPV[who,"Lower"],upper=object$CI.PPV[who,"Upper"]),NPV=round(myFactor*prNPV,myDigits),"(CI.NPV)"=MgFormCi(lower=object$CI.NPV[who,"Lower"],upper=object$CI.NPV[who,"Upper"])))
##     LRplus <- prSens/(1-prSpec)
##     LRminus <- (1-prSens)/prSpec
##     if (length(breaks)==NROW(out)-1)
##       rownames(out) <- c(round(breaks,digits),"--")
##     else
##       rownames(out) <- round(breaks,digits)
##     cat("\n\nSens=Sensitivity\nSpec=Specificity\nPPV=positive predictive value\nNPV=negative predictive value\n\n")
##     print(out,quote=FALSE)
##     cat("\n\nLR+=positive likelihood ratio\nLR-=negative likelihood ratio\n\n")
##     print(cbind(rownames(out),"LR+"=round(LRplus,digits),"LR-"=round(LRminus,digits)),quote=FALSE)
##   }
##   else{
##     out <- do.call("cbind",list(Sens=object$Roc$Sensitivity[who],Spec=object$Roc$Specificity[who],PPV=object$Roc$PPV[who],NPV=object$Roc$NPV[who]))
##     rownames(out) <- round(breaks,digits)
##     cat("\n\nSens=Sensitivity\nSpec=Specificity\nPPV=positive predictive value\nNPV=negative predictive value\n\n")
##     print(100*out,digits=digits)
##   }
## }
