rm(list=ls())

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(OCNet)
library(spam)
library(tictoc)
library(optimParallel)
source("eDITH_functions.R")

load("OCN_eDNA_SI.rda") # this loads the OCn with Q0=2.5 m3s-1. Change to "OCN_eDNA.rda" for the version with Q0=10m3s-1
OCN$RN$downstreamPath <- NULL # save memory

tau_sol <- 4*3600 # change to 1*3600 for the version with tau=1 h
filestring <- "results_tau4_Q25.rda" # change name accordingly if a different version is run

modelType <- c("me")
strategyType <- c("u0","u20", "u50", "u80", "u100")
samplingIntensity <- c("s10","s25","s50") 
distributionType <- c("R1","R2","R3","R4","R5","A1","A2","A3","A4","A5")

Nodes <- setdiff(1:OCN$AG$nNodes,OCN$AG$outlet)
NodesDownstream <- intersect(which(OCN$AG$A > median(OCN$AG$A[Nodes])), Nodes)
NodesUpstream <- intersect(which(OCN$AG$A <= median(OCN$AG$A[Nodes])), Nodes)

# define hotspots ####
availableNodes <- setdiff(1:OCN$AG$nNodes, OCN$AG$outlet)
hotspotList <- vector("list",0)
set.seed(1); hotspotList$R1 <- sample(availableNodes, 5)
set.seed(2); hotspotList$R2 <- sample(availableNodes, 5)
set.seed(3); hotspotList$R3 <- sample(availableNodes, 5)
set.seed(4); hotspotList$R4 <- sample(availableNodes, 5)
set.seed(5); hotspotList$R5 <- sample(availableNodes, 5)
set.seed(11); hotspotList$A1 <- sample(availableNodes, 50)
set.seed(12); hotspotList$A2 <- sample(availableNodes, 50)
set.seed(13); hotspotList$A3 <- sample(availableNodes, 50)
set.seed(14); hotspotList$A4 <- sample(availableNodes, 50)
set.seed(15); hotspotList$A5 <- sample(availableNodes, 50)


# pick sampling sites ####
samplingSitesList <- vector("list",length(samplingIntensity))
names(samplingSitesList) <- samplingIntensity
k <- 0
for (sampling in samplingIntensity){
  nSampling <- as.numeric(gsub("[^0-9.]", "",  sampling))*2
  for (strategy in strategyType){
    samplingSitesList[[sampling]][[strategy]] <- matrix(0,ncol=nSampling,nrow=10)
    
    weightUpstream <- as.numeric(gsub("[^0-9.]", "",  strategy))/100
    weightDownstream <- 1 - weightUpstream
    
    for (indSim in 1:10){
      k <- k+1
      set.seed(k) # seed changing everytime
      samplingSitesList[[sampling]][[strategy]][indSim,] <- c( sample(NodesDownstream,round(weightDownstream*nSampling)), 
                                                               sample(NodesUpstream,round(weightUpstream*nSampling))) 
      
    }
  }
}

# perform simulations ####
results_all <- vector("list",length(modelType))
names(results_all) <- modelType

# open parallel clusters
cl <- makeCluster(spec=detectCores(), outfile="")
setDefaultCluster(cl=cl)
clusterExport(cl, c("eval_conc","eDITH"))
clusterEvalQ(cl, library("OCNet"))

set.seed(200) # this seed controls initialParamValues 

for (model in modelType){
  for (strategy in strategyType){
    for (sampling in samplingIntensity){
      for (distrib in distributionType){ 
        for (indSim in 1:10){
          
          p_sol <- numeric(OCN$AG$nNodes)
          
          hotspots <- hotspotList[[distrib]]
          halfHotspots <- vector("numeric",0)
          for (ih in hotspots){
            halfHotspots <- c(halfHotspots, OCN$AG$downNode[ih]) 
            halfHotspots <- c(halfHotspots, which(OCN$AG$downNode==ih))
            p_sol[ih] <- p_sol[ih] + 1
          }
          for (ihh in halfHotspots){
            p_sol[ihh] <- p_sol[ihh] + 0.5
          }
          p_sol <- p_sol/sum(p_sol)*sum(OCN$AG$leng * OCN$AG$width) / (OCN$AG$leng * OCN$AG$width)
          p_sol[OCN$AG$outlet] <- 0
          
          C_sol <- eval_conc(OCN, tau_sol, p_sol, "AG", normalize=TRUE)
          
          # add measurement error
          if (model == "me"){
            C_obs <- C_sol
            C_obs[runif(length(C_sol)) < exp(-C_sol/1)] <- 0 # non detection probability 
            C_obs <- exp(rnorm(length(C_obs),log(C_obs),0.5)) 
          } else {
            C_obs <- C_sol
          }
          
          # choose sampling sites
          samplingSites <- samplingSitesList[[sampling]][[strategy]][indSim,]
          
          cat("\n")
          cat(sprintf("sim: %d    -   model: %s    -   strategy: %s   -   intensity: %s   -   distrib: %s   -   %11s \n",
                      indSim,model,strategy,sampling,distrib,format(Sys.time(),"%b%d %H:%M")))
          
          # fit model
          if (model == "tf"){
            initialParamValues <- runif(OCN$AG$nNodes,0,5) 
            optList <- optimParallel(initialParamValues, eval_loglik_taufixed, 
                                     lower=0,upper=250,
                                     tau_sol=tau_sol, samplingSites=samplingSites, ConcAG=C_obs, OCN=OCN,
                                     control=list(fnscale=-1, maxit=200,trace=1, REPORT=10))  
          } else {
            initialParamValues <- c(runif(1,0,7), runif(OCN$AG$nNodes,0,5)) 
            optList <- optimParallel(initialParamValues, eval_loglik, 
                                     lower=0,upper=250,
                                     samplingSites=samplingSites, ConcAG=C_obs, OCN=OCN,
                                     control=list(fnscale=-1, maxit=200,trace=1, REPORT=10))
          }
          
          # save results
          results_all[[model]][[strategy]][[sampling]][[distrib]][[as.character(indSim)]][["par"]] <- optList$par
          results_all[[model]][[strategy]][[sampling]][[distrib]][[as.character(indSim)]][["samplingSites"]] <- samplingSites
          results_all[[model]][[strategy]][[sampling]][[distrib]][[as.character(indSim)]][["C_sol"]] <- C_sol
          results_all[[model]][[strategy]][[sampling]][[distrib]][[as.character(indSim)]][["C_obs"]] <- C_obs
          results_all[[model]][[strategy]][[sampling]][[distrib]][[as.character(indSim)]][["p_sol"]] <- p_sol
          results_all[[model]][[strategy]][[sampling]][[distrib]][[as.character(indSim)]][["value"]] <- optList$value
          results_all[[model]][[strategy]][[sampling]][[distrib]][[as.character(indSim)]][["convergence"]] <- optList$convergence
          results_all[[model]][[strategy]][[sampling]][[distrib]][[as.character(indSim)]][["value_sol"]] <- eval_loglik(c(tau_sol/3600 - 1, p_sol), samplingSites, C_sol, OCN)
          if (model == "tf"){
            results_all[[model]][[strategy]][[sampling]][[distrib]][[as.character(indSim)]][["C"]] <- eval_conc(OCN, tau_sol, optList$par, "AG", normalize=TRUE)  
          } else {
            results_all[[model]][[strategy]][[sampling]][[distrib]][[as.character(indSim)]][["C"]] <- eval_conc(OCN, (optList$par[1] + 1)*3600, optList$par[-1], "AG", normalize=TRUE)
          }
          results_all[[model]][[strategy]][[sampling]][[distrib]][[as.character(indSim)]][["initialParamValues"]] <- initialParamValues
          
          save(results_all,file=filestring)
        }
      }
    }
  }
}
setDefaultCluster(cl=NULL); stopCluster(cl)
