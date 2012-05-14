getResponse <- function(formula,refLevel,data,verbose){
  ## case 1: continuous 
  ## case 2: binary
  ## case 3: ordinal
  ## case 4: survival (Surv, Hist) 
  ## case 5: competing risks (Hist)
  ## case 6: updating (?)
  responseVariables <- all.vars(formula)
  if (length(responseVariables)==1){
    ## case 1,2,3
    m <- model.frame(formula=formula,data=data,na.action=na.fail)
    response <- model.response(m)
    if (is.factor(response) || length(unique(response))==2){
      response <- factor(response)
      if (length(levels(response))==2) {
        if (missing(refLevel)) refLevel <- levels(response)[2]
        ## response <- as.numeric(response==refLevel)
        attr(response,"refLevel") <- refLevel
        attr(response,"model") <- "binary"
      }
      else{
        attr(response,"model") <- "multi.level"
        stop("Methods for factors with more than two classes are not\n (not yet) implemented.")
      }
    }
    else{
      attr(response,"model") <- "continuous"
    }
  }
  else{
    responseNames <- all.names(formula)
    if (responseNames[2] %in% c("Surv","Hist")){
      ## case 4,5
      ## replace Surv by Hist to use the attributes of
      ## the Hist object
      if (responseNames[2]=="Surv"){
        formula[[2]][[1]] <- as.name("Hist")
      }
      m <- model.frame(formula=formula,data=data,na.action=na.fail)
      response <- model.response(m)
      ## if (attr(response,"model")=="competing.risks"){
      ## if (missing(refLevel)) refLevel <- attr(response,"states")[1]
      ## attr(response,"cause") <- refLevel
      ## }
    }
    else{
      stop("Cannot assign response type.")
    }
  }
  response
}
