#function to get subjective predicted probabilities of choice for each choice option
predictIndProb <-function(X, DrawCoeff, deData){
  eXb <- as.numeric(exp(X %*% DrawCoeff))  #gettin exp 
  #SeXb <- tapply(eXb, deData[,"chid"], sum)
  SeXb <-rowsum(eXb, deData[,"chid"])
  P <- eXb / SeXb[deData[,"chid"]]
  return(P)
}

calcK <-function(iterData, i){
  p_index <- which(colnames(iterData)==paste("l_", i, sep=""))
  kbb <- as.numeric(iterData[,'cVar'] * iterData[,p_index])
  kbb <- ifelse(kbb > 0, kbb, NA)
  iterData[, paste("k", i, sep="")] <- kbb
  WM[, paste("k", i, sep="")] <<- kbb
  kbbb <- prod(iterData[,paste("k", i, sep="")], na.rm = T)
  return(kbbb)  
}


# Pnobs ->obs foreach people
# k -> number of normals
# options
#
createWM <- function(Pnobs, k, options, probs){
  #--------------------------------
  # Create working matrix
  #--------------------------------
  WM <- matrix(, nrow=Pnobs)
  myvector<-c(NA)
  length(myvector)<-length(Pnobs)
  
  cols <- c("l_", "prob_", "k")
  for(c in cols){
    for(n in 1:k){
      WM <- cbind(WM, colName= myvector)
      colnames(WM)[length(colnames(WM))]<- paste(c, n, sep="")
    }
  }
  WM <-WM[,2:length(WM[1,])] #delete first column created for memory pre-allocation
  #create grouping column
  WM<-cbind(WM, "chid"=rep(1:(nrow(WM)/options), each=options))
  WM<-cbind(WM, "cVar"=myvector)
  #c_index <- which(colnames(WM)=="cVar")
  #set initial probability & logit
  for(i in 1:k){
    WM[, paste("prob_", i, sep="")] <- probs[i]
  }
  return(WM)
}

createResM <- function(Pnobs, k, R, nobs){
  myvector<-c(NA)
  length(myvector)<-length(Pnobs)
  
  resData <- matrix(, nrow=nobs*R)
  resData<-cbind(resData,'subj'=rep(1:nobs, each=R))
  resData<-cbind(resData,'draw'=rep(1:R, times=nobs))
  resData <- cbind(resData, "den" = NA)
  resData <-resData[,2:length(resData[1,])] #delete first column created for memory pre-allocation
  cols <- c("K", "h")
  for(c in cols){
    for(n in 1:k){
      resData <- cbind(resData, colName= myvector)
      colnames(resData)[length(colnames(resData))]<- paste(c, n, sep="")
    }
  }
  
  return(resData)
}

getExecTotalTime <- function(totalTime){
  #get total execution time
  totalTimeElapsed <- proc.time() - totalTime
  myTime <- round(totalTimeElapsed["elapsed"], 2)
  minutes <-as.integer(myTime/60)
  secs <- round(myTime - minutes*60, 0)
  
  cat(sprintf("\nTotal execution m: %s    s: %s \n\n\n", minutes, secs))
}