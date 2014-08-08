setwd("~/Google Drive/Universita/R/scripts/biomasse_ce/data_analysis/EM")
library('randtoolbox') #for halton draws
#get execution time
totalTime <- proc.time()

#load functions
source('EM_mixedNormals_functions.R')

#--------------------------------
# DATA
#--------------------------------
#get data
datasetFinal <- read.csv("../datasetBiomasse.csv", header = TRUE, sep = ",")
datasetFinal <-datasetFinal[datasetFinal$protest==0,]
attach(datasetFinal)
EMData <- cbind(subj, set_id, altern, altern, size_s,size_b,dist3s,dist10s,dist3b,dist10b,cert, cost,mychoice)
detach(datasetFinal)
EMData<-as.data.frame(EMData)

#keep only nobs respondents of full dataset
EMData<-EMData[1:(nobs*options*cset),]

#--------------------------------
# DEFINE NEEDED values
#--------------------------------
nIter = 2 #set number of iterations
R = 10 #number of draws per person
options=4 # num of choice alternatives per choice set
cset <- 6 #number of choice sets per person

nobs=5 #total number of persons (subjects) in dataset
dim=9; #dimension of simulated data  (betas)
k=2  #number of "true" components of mixed normals


#choice data vars
id<-"subj"        # variable that indentifies each respondent
cVar<-"mychoice"  # choice variable
optVar<-"option"  # var that identifies the options in each choice set
#--------------------------------
# END INPUTS
#-------------------------------



Pnobs <- options*cset #total number of observations per person
#create grouping column
EMData<-cbind(EMData, "chid"=rep(1:(nrow(EMData)/options), each=options))

#--------------------------------
# Initial probability of belonging to normal class 
#--------------------------------
Sc <- vector()
for(i in 1:k){
  Sc[i] <- 1/k 
}


#--------------------------------
# Create matrix for storing results
#--------------------------------
resData <- createResM(Pnobs, k, R, nobs)

#---------------
# MAKE random DRAWS for k-th normal for N-th subject
#---------------
ParMatrix <- matrix(, nrow=R, ncol=k*dim)

#W <-matrix(0, nrow=9, ncol=9)
#diag(W)<-1

diffLL <- 1
storedParMeans<-matrix(nrow=0, ncol=k*dim)
parMeans <-rep(0, times=k*dim)
parCov_trans<-list()#for storing iteration vals
parCov <-rep(1, times=k*dim)
parCov_trans[[1]]<-parCov

pCov<-list()
#storedParMeans<-rbind(storedParMeans, parMeans)

Hseq <- halton(R + 10, dim=18, normal=TRUE) #halton draw
#dischard first 10 draws (see Train)
Hseq<-Hseq[11:nrow(Hseq),]
iterLL <-c(NA)
iterLL<-length(nIter)
Iter <- 0
#set convergence cotroller
converged <-FALSE
while(Iter <= nIter && converged!=TRUE){
  Rprof()
  #get starting time
  ptm <- proc.time()
  #--------------------------------
  # Create matrix for calulations during iterations
  #--------------------------------
  WM <- createWM(Pnobs, k, options, Sc)
  
  totalPar <-matrix(nrow=R*nobs, ncol=k*dim)
  #loop counter
  Rcounter = 0
  
  #calculate
  for(n in 1:nobs){
    parStart <-((n-1)*R+1)
    parEnd<-n*R
    #make drwas for subject
    for (i in 1:(k*dim)){
      #Bi <-parMeans[i] + rnorm(R, 0, 1)
      #Bi <-rnorm(R, parMeans[i], sqrt(parCov_trans[[Iter+1]][i])
        
      Bi<-parMeans[i] + parCov_trans[[Iter+1]][i]*Hseq[,i]
      ParMatrix[,i] <-Bi
    }
    totalPar[parStart:parEnd,] <-ParMatrix
    subjStart <-((n-1)*Pnobs+1)
    subjEnd<-n*Pnobs
    X <- as.matrix(EMData[subjStart:subjEnd, 1:dim]) #matrix of X for subject n @todo --> take out of loop for speed
    #foreach draw of params
    for(Ri in 1:R){
      Rcounter = Rcounter + 1
      for(Kn in 1:k){
        mStart <-((Kn-1)*dim+1)
        mEnd<-Kn*dim
        draw <- ParMatrix[Ri, mStart:mEnd]
        
        WM[, paste("l_", Kn, sep="")] <- predictIndProb(X, draw, WM)
        WM[, 'cVar']<-EMData[subjStart:subjEnd,'mychoice']
        #WM[Ri, paste("kbb", Kn, sep="")]<-calcKbbb(WM, Kn)
        resData[Rcounter,paste("K", Kn, sep="")] <-calcK(WM, Kn)
      }
    }
    #Pn<-
  }
  resData[,"den"]<-(resData[,4:(4+k-1)]%*%Sc)/R
  for(Kn in 1:k){
    resData[,paste("h", Kn, sep="")] <-Sc[Kn]*resData[,paste("K", Kn, sep="")]/resData[,"den"]
    # UPDATE SHARES Sn
    Sc[Kn] <-sum(resData[,paste("h", Kn, sep="")])
  }
  
  Stot <-sum(Sc)
  Sc<-Sc/Stot
  
  sumLL <- round(sum(log(resData[,"den"]), na.rm = T), 4)
  iterLL[Iter]<-sumLL
    
    #----------------------------------------
    # calculate new par means and covariances
    #----------------------------------------
  for(Kn in 1:k){
      mStart <-((Kn-1)*dim+1)
      mEnd<-Kn*dim
      parMeans[mStart:mEnd] <- apply(totalPar[, mStart:mEnd] * resData[,paste("h", Kn, sep="")],2,sum)/sum(resData[,paste("h", Kn, sep="")])
      Bdiff <-totalPar[, mStart:mEnd] - parMeans[mStart:mEnd]
      if(Kn==1){
        parCov1<-cov(Bdiff* resData[,paste("h", Kn, sep="")]/sum(resData[,paste("h", Kn, sep="")]))
      }else{
        parCov2<-cov(Bdiff* resData[,paste("h", Kn, sep="")]/sum(resData[,paste("h", Kn, sep="")]))
      }
      #parCov<-apply((Bdiff*I(Bdiff)) * resData[,paste("h", Kn, sep="")],2,sum)/sum(resData[,paste("h", Kn, sep="")])
      #pippo<- cov.wt(totalPar[, mStart:mEnd], resData[,paste("h", Kn, sep="")], method = "ML", cor = TRUE)
      #pCov[[Kn]]<-diag(chol(pippo$cov))
      #parCov_trans[Kn]<- apply((Bdiff*I(Bdiff)) * resData[,paste("h", Kn, sep="")],2,sum)/sum(resData[,paste("h", Kn, sep="")])
    
    }
    parCov_trans[[Iter+2]] <-c(parCov1, parCov2)
    #parCov_trans[[Iter+2]] <-c(pCov[[1]], pCov[[2]])
    #storedParMeans<-rbind(storedParMeans, parMeans)
  
  Rprof(NULL)
  timeElapsed <- proc.time() - ptm
  myTime <- round(timeElapsed["elapsed"], 2)
  
  cat(sprintf("Iteration %s log likelihood = %s    time:%s sec\n", Iter, sumLL, myTime))
  
  #----------------------------------------------
  # exit from iterations if log likelihood is the same in last iteration
  #----------------------------------------------
  
  if(Iter >=2){
    comparison <- Iter - 1
  }
  
  if(Iter > 2 && {iterLL[Iter]-iterLL[comparison]} > 0 && {iterLL[Iter]-iterLL[comparison]}< 0.00005){
    
    cat("------------------------\n")
    cat(sprintf("CONVERGENCE REACHED @ iteration  %s\n", Iter))
    cat("------------------------\n")
    converged <-TRUE
    
  }
  else if(converged==FALSE && Iter==nIter){
    
      cat("------------------------\n")
      cat("CONVERGENCE NOT REACHED: consider increasing the numbers of iterations \n")
      cat("------------------------\n")
    #increase counter
    
    
  }
  
  plot(iterLL, type="l")
  
  Iter <- Iter + 1 #increase iteration counter
}

getExecTotalTime(totalTime)
prof = summaryRprof()
prof$by.self[1:5,]