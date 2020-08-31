rm(list=ls())

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

draw_figs <- 0 # set to 1 to draw figs

library("spam")
library("OCNet")

# load data ####

source("eDITH_functions.R")
load("OCN_eDNA.rda")
OCN$RN$downstreamPath <- NULL

simulationType <- c("me") # 
strategyType <- c("u0","u20","u50","u80","u100") 
samplingIntensity <- c("s10","s25","s50") #
distributionType <- c("R","A")
distributionNo <- c(1:5)
nRepeat <- 10
nTotSim <- length(simulationType)*length(strategyType)*length(samplingIntensity)*length(distributionType)*length(distributionNo)*nRepeat

outlet <- OCN$AG$outlet

Nodes <- setdiff(1:OCN$AG$nNodes,OCN$AG$outlet)
NodesDownstream <- intersect(which(OCN$AG$A > median(OCN$AG$A[Nodes])), Nodes)
NodesUpstream <- intersect(which(OCN$AG$A <= median(OCN$AG$A[Nodes])), Nodes)
nDwn <- NodesDownstream
nDwn <- nDwn[nDwn==64]


simulation <- vector("character",nTotSim*4)
strategy <- distribution_type <- distribution_no <- simulation
intensity <- vector("numeric",nTotSim)
loose <- strict <- p_a <- mab_p <- mab_C <- tau <- tau_sol <- Q <- intensity
looseDown <- looseUp <- strictDown <- strictUp  <- p_aDown <- p_aUp <- intensity


# read results for tau=4 h, Q0=2.5 m3/s ####
load("results_tau4_Q25.rda")
k <- 0
for (indSim in 1:length(simulationType)){ #
  sim <- simulationType[indSim]
  for (indStrategy in 1:length(strategyType)){
    strat <- strategyType[indStrategy]
    for (indSampling in 1:length(samplingIntensity)){
      sampl <- samplingIntensity[indSampling]
      for (indDistrib in 1:length(distributionType)){
        distr_t <- distributionType[indDistrib]
        for (distr_n in 1:length(distributionNo)){
          distr <- paste(distr_t,as.character(distr_n),sep="")
          for (indRep in 1:nRepeat){
            k <- k + 1
            res <- results_all[[sim]][[strat]][[as.character(sampl)]][[distr]][[as.character(indRep)]]
            
            if (!is.null(res)){
              param <- res$par
              if (sim != "tf"){
                p <- param[-1] # remove tau
              } else {p <- param}
              p[outlet] <- NaN # remove outlet
              p_sol <- res$p_sol
              p_sol[outlet] <- NaN
              C_sol <- res$C_sol
              C <- res$C_obs
              simulation[k] <- sim
              strategy[k] <- strat
              intensity[k] <- sampl
              distribution_type[k] <- distr_t
              distribution_no[k] <- distr_n
              
              loose[k] <- sum( (p_sol <= 0.2 & p <= 0.2) | (p >= 0.5*p_sol & p <= 2*p_sol), na.rm=T) / 200
              strict[k] <- sum( (p_sol <= 0.25 & p <= 0.25) | (p >= 0.75*p_sol & p <= 4/3*p_sol), na.rm=T) / 200
              p_a[k] <- sum( (p_sol <= 0.25 & p <= 0.25) | (p_sol > 0.2 & p > 0.2), na.rm=T) / 200
              mab_p[k] <- mean(abs(p_sol - p), na.rm=T)
              mab_C[k] <- mean(abs(C_sol - C))
              
              looseDown[k] <- sum( (p_sol[NodesDownstream] <= 0.2 & p[NodesDownstream] <= 0.2) |
                                     (p[NodesDownstream] >= 0.5*p_sol[NodesDownstream] & p[NodesDownstream] <= 2*p_sol[NodesDownstream]) )/length(NodesDownstream)
              looseUp[k] <- sum( (p_sol[NodesUpstream] <= 0.2 & p[NodesUpstream] <= 0.2) |
                                   (p[NodesUpstream] >= 0.5*p_sol[NodesUpstream] & p[NodesUpstream] <= 2*p_sol[NodesUpstream]) )/length(NodesDownstream)
              strictDown[k] <- sum( (p_sol[NodesDownstream] <= 0.2 & p[NodesDownstream] <= 0.2) |
                                      (p[NodesDownstream] >= 0.75*p_sol[NodesDownstream] & p[NodesDownstream] <= 4/3*p_sol[NodesDownstream]) )/length(NodesDownstream)
              strictUp[k] <- sum( (p_sol[NodesUpstream] <= 0.2 & p[NodesUpstream] <= 0.2) |
                                    (p[NodesUpstream] >= 0.75*p_sol[NodesUpstream] & p[NodesUpstream] <= 4/3*p_sol[NodesUpstream]) )/length(NodesDownstream)
              p_aDown[k] <- sum( (p_sol[NodesDownstream] <= 0.2 & p[NodesDownstream] <= 0.2) |
                                   (p[NodesDownstream] >= 0.5*p_sol[NodesDownstream] & p[NodesDownstream] <= 2*p_sol[NodesDownstream]) )/length(NodesDownstream)
              p_aUp[k] <- sum( (p_sol[NodesUpstream] <= 0.2 & p[NodesUpstream] <= 0.2) |
                                 (p_sol[NodesUpstream] > 0.2 & p[NodesUpstream] > 0.2) )/length(NodesDownstream)
              tau_sol[k] <- 4
                Q[k] <- 25
                
                if (sim=="tf"){
                  tau[k] <- 4
                } else {tau[k] <- 1 + param[1]}
              
            }
          }
        }
      }
    }
  }
}

# read results for tau=1 h, Q0=10 m3/s ####
load("results_tau1_Q10.rda")
for (indSim in 1:length(simulationType)){ #
  sim <- simulationType[indSim]
  for (indStrategy in 1:length(strategyType)){
    strat <- strategyType[indStrategy]
    for (indSampling in 1:length(samplingIntensity)){
      sampl <- samplingIntensity[indSampling]
      for (indDistrib in 1:length(distributionType)){
        distr_t <- distributionType[indDistrib]
        for (distr_n in 1:length(distributionNo)){
          distr <- paste(distr_t,as.character(distr_n),sep="")
          for (indRep in 1:nRepeat){
            k <- k + 1
            res <- results_all[[sim]][[strat]][[as.character(sampl)]][[distr]][[as.character(indRep)]]
            
            if (!is.null(res)){
              param <- res$par
              if (sim != "tf"){
                p <- param[-1] # remove tau
              } else {p <- param}
              p[outlet] <- NaN # remove outlet
              p_sol <- res$p_sol
              p_sol[outlet] <- NaN
              C_sol <- res$C_sol
              C <- res$C_obs
              simulation[k] <- sim
              strategy[k] <- strat
              intensity[k] <- sampl
              distribution_type[k] <- distr_t
              distribution_no[k] <- distr_n
              
              loose[k] <- sum( (p_sol <= 0.2 & p <= 0.2) | (p >= 0.5*p_sol & p <= 2*p_sol), na.rm=T) / 200
              strict[k] <- sum( (p_sol <= 0.25 & p <= 0.25) | (p >= 0.75*p_sol & p <= 4/3*p_sol), na.rm=T) / 200
              p_a[k] <- sum( (p_sol <= 0.25 & p <= 0.25) | (p_sol > 0.2 & p > 0.2), na.rm=T) / 200
              mab_p[k] <- mean(abs(p_sol - p), na.rm=T)
              mab_C[k] <- mean(abs(C_sol - C))
              
              looseDown[k] <- sum( (p_sol[NodesDownstream] <= 0.2 & p[NodesDownstream] <= 0.2) |
                                     (p[NodesDownstream] >= 0.5*p_sol[NodesDownstream] & p[NodesDownstream] <= 2*p_sol[NodesDownstream]) )/length(NodesDownstream)
              looseUp[k] <- sum( (p_sol[NodesUpstream] <= 0.2 & p[NodesUpstream] <= 0.2) |
                                   (p[NodesUpstream] >= 0.5*p_sol[NodesUpstream] & p[NodesUpstream] <= 2*p_sol[NodesUpstream]) )/length(NodesDownstream)
              strictDown[k] <- sum( (p_sol[NodesDownstream] <= 0.2 & p[NodesDownstream] <= 0.2) |
                                      (p[NodesDownstream] >= 0.75*p_sol[NodesDownstream] & p[NodesDownstream] <= 4/3*p_sol[NodesDownstream]) )/length(NodesDownstream)
              strictUp[k] <- sum( (p_sol[NodesUpstream] <= 0.2 & p[NodesUpstream] <= 0.2) |
                                    (p[NodesUpstream] >= 0.75*p_sol[NodesUpstream] & p[NodesUpstream] <= 4/3*p_sol[NodesUpstream]) )/length(NodesDownstream)
              p_aDown[k] <- sum( (p_sol[NodesDownstream] <= 0.2 & p[NodesDownstream] <= 0.2) |
                                   (p[NodesDownstream] >= 0.5*p_sol[NodesDownstream] & p[NodesDownstream] <= 2*p_sol[NodesDownstream]) )/length(NodesDownstream)
              p_aUp[k] <- sum( (p_sol[NodesUpstream] <= 0.2 & p[NodesUpstream] <= 0.2) |
                                 (p_sol[NodesUpstream] > 0.2 & p[NodesUpstream] > 0.2) )/length(NodesDownstream)
              tau_sol[k] <- 1
                Q[k] <- 10
                
                if (sim=="tf"){
                  tau[k] <- 4
                } else {tau[k] <- 1 + param[1]}
              
            }
          }
        }
      }
    }
  }
}

# read results for tau=1 h, Q0=2.5 m3/s ####
load("results_tau1_Q25.rda")
for (indSim in 1:length(simulationType)){ #
  sim <- simulationType[indSim]
  for (indStrategy in 1:length(strategyType)){
    strat <- strategyType[indStrategy]
    for (indSampling in 1:length(samplingIntensity)){
      sampl <- samplingIntensity[indSampling]
      for (indDistrib in 1:length(distributionType)){
        distr_t <- distributionType[indDistrib]
        for (distr_n in 1:length(distributionNo)){
          distr <- paste(distr_t,as.character(distr_n),sep="")
          for (indRep in 1:nRepeat){
            k <- k + 1
            res <- results_all[[sim]][[strat]][[as.character(sampl)]][[distr]][[as.character(indRep)]]
            
            if (!is.null(res)){
              param <- res$par
              if (sim != "tf"){
                p <- param[-1] # remove tau
              } else {p <- param}
              p[outlet] <- NaN # remove outlet
              p_sol <- res$p_sol
              p_sol[outlet] <- NaN
              C_sol <- res$C_sol
              C <- res$C_obs
              simulation[k] <- sim
              strategy[k] <- strat
              intensity[k] <- sampl
              distribution_type[k] <- distr_t
              distribution_no[k] <- distr_n
              
              loose[k] <- sum( (p_sol <= 0.2 & p <= 0.2) | (p >= 0.5*p_sol & p <= 2*p_sol), na.rm=T) / 200
              strict[k] <- sum( (p_sol <= 0.25 & p <= 0.25) | (p >= 0.75*p_sol & p <= 4/3*p_sol), na.rm=T) / 200
              
              p_a[k] <- sum( (p_sol <= 0.25 & p <= 0.25) | (p_sol > 0.2 & p > 0.2), na.rm=T) / 200
              results_all[[sim]][[strat]][[as.character(sampl)]][[distr]][[as.character(indRep)]]$p_a <- p_a[k]
              
              mab_p[k] <- mean(abs(p_sol - p), na.rm=T)
              mab_C[k] <- mean(abs(C_sol - C))
              
              looseDown[k] <- sum( (p_sol[NodesDownstream] <= 0.2 & p[NodesDownstream] <= 0.2) |
                                     (p[NodesDownstream] >= 0.5*p_sol[NodesDownstream] & p[NodesDownstream] <= 2*p_sol[NodesDownstream]) )/length(NodesDownstream)
              looseUp[k] <- sum( (p_sol[NodesUpstream] <= 0.2 & p[NodesUpstream] <= 0.2) |
                                   (p[NodesUpstream] >= 0.5*p_sol[NodesUpstream] & p[NodesUpstream] <= 2*p_sol[NodesUpstream]) )/length(NodesDownstream)
              strictDown[k] <- sum( (p_sol[NodesDownstream] <= 0.2 & p[NodesDownstream] <= 0.2) |
                                      (p[NodesDownstream] >= 0.75*p_sol[NodesDownstream] & p[NodesDownstream] <= 4/3*p_sol[NodesDownstream]) )/length(NodesDownstream)
              strictUp[k] <- sum( (p_sol[NodesUpstream] <= 0.2 & p[NodesUpstream] <= 0.2) |
                                    (p[NodesUpstream] >= 0.75*p_sol[NodesUpstream] & p[NodesUpstream] <= 4/3*p_sol[NodesUpstream]) )/length(NodesDownstream)
              p_aDown[k] <- sum( (p_sol[NodesDownstream] <= 0.2 & p[NodesDownstream] <= 0.2) |
                                   (p[NodesDownstream] >= 0.5*p_sol[NodesDownstream] & p[NodesDownstream] <= 2*p_sol[NodesDownstream]) )/length(NodesDownstream)
              p_aUp[k] <- sum( (p_sol[NodesUpstream] <= 0.2 & p[NodesUpstream] <= 0.2) |
                                 (p_sol[NodesUpstream] > 0.2 & p[NodesUpstream] > 0.2) )/length(NodesDownstream)
              tau_sol[k] <- 1 
                Q[k] <- 25
                
                if (sim=="tf"){
                  tau[k] <- 4
                } else {tau[k] <- 1 + param[1]}
              
            }
          }
        }
      }
    }
  }
}

# read results for tau=4 h, Q0=10 m3/s ####
load("results_all.rda")
for (indSim in 1:length(simulationType)){ #
  sim <- simulationType[indSim]
  for (indStrategy in 1:length(strategyType)){
    strat <- strategyType[indStrategy]
    for (indSampling in 1:length(samplingIntensity)){
      sampl <- samplingIntensity[indSampling]
      for (indDistrib in 1:length(distributionType)){
        distr_t <- distributionType[indDistrib]
        for (distr_n in 1:length(distributionNo)){
          distr <- paste(distr_t,as.character(distr_n),sep="")
          for (indRep in 1:nRepeat){
            k <- k + 1
            res <- results_all[[sim]][[strat]][[as.character(sampl)]][[distr]][[as.character(indRep)]]
            
            if (!is.null(res)){
              param <- res$par
              if (sim != "tf"){
                p <- param[-1] # remove tau
              } else {p <- param}
              p[outlet] <- NaN # remove outlet
              p_sol <- res$p_sol
              p_sol[outlet] <- NaN
              C_sol <- res$C_sol
              C <- res$C_obs
              simulation[k] <- sim
              strategy[k] <- strat
              intensity[k] <- sampl
              distribution_type[k] <- distr_t
              distribution_no[k] <- distr_n
              
              loose[k] <- sum( (p_sol <= 0.2 & p <= 0.2) | (p >= 0.5*p_sol & p <= 2*p_sol), na.rm=T) / 200
              strict[k] <- sum( (p_sol <= 0.25 & p <= 0.25) | (p >= 0.75*p_sol & p <= 4/3*p_sol), na.rm=T) / 200
              results_all[[sim]][[strat]][[as.character(sampl)]][[distr]][[as.character(indRep)]]$strict <- strict[k]
              
              p_a[k] <- sum( (p_sol <= 0.25 & p <= 0.25) | (p_sol > 0.2 & p > 0.2), na.rm=T) / 200
              results_all[[sim]][[strat]][[as.character(sampl)]][[distr]][[as.character(indRep)]]$p_a <- p_a[k]
              
              mab_p[k] <- mean(abs(p_sol - p), na.rm=T)
              mab_C[k] <- mean(abs(C_sol - C))
              
              looseDown[k] <- sum( (p_sol[NodesDownstream] <= 0.2 & p[NodesDownstream] <= 0.2) |
                                     (p[NodesDownstream] >= 0.5*p_sol[NodesDownstream] & p[NodesDownstream] <= 2*p_sol[NodesDownstream]) )/length(NodesDownstream)
              looseUp[k] <- sum( (p_sol[NodesUpstream] <= 0.2 & p[NodesUpstream] <= 0.2) |
                                   (p[NodesUpstream] >= 0.5*p_sol[NodesUpstream] & p[NodesUpstream] <= 2*p_sol[NodesUpstream]) )/length(NodesDownstream)
              strictDown[k] <- sum( (p_sol[NodesDownstream] <= 0.2 & p[NodesDownstream] <= 0.2) |
                                      (p[NodesDownstream] >= 0.75*p_sol[NodesDownstream] & p[NodesDownstream] <= 4/3*p_sol[NodesDownstream]) )/length(NodesDownstream)
              strictUp[k] <- sum( (p_sol[NodesUpstream] <= 0.2 & p[NodesUpstream] <= 0.2) |
                                    (p[NodesUpstream] >= 0.75*p_sol[NodesUpstream] & p[NodesUpstream] <= 4/3*p_sol[NodesUpstream]) )/length(NodesDownstream)
              p_aDown[k] <- sum( (p_sol[NodesDownstream] <= 0.2 & p[NodesDownstream] <= 0.2) |
                                   (p[NodesDownstream] >= 0.5*p_sol[NodesDownstream] & p[NodesDownstream] <= 2*p_sol[NodesDownstream]) )/length(NodesDownstream)
              p_aUp[k] <- sum( (p_sol[NodesUpstream] <= 0.2 & p[NodesUpstream] <= 0.2) |
                                 (p_sol[NodesUpstream] > 0.2 & p[NodesUpstream] > 0.2) )/length(NodesDownstream)
              tau_sol[k] <- 4
              Q[k] <- 10
              
              if (sim=="tf"){
                tau[k] <- 4
              } else {tau[k] <- 1 + param[1]}
              
            }
          }
        }
      }
    }
  }
}

# create dataframe
data <- data.frame(simulation,strategy,intensity,distribution_type,distribution_no,tau_sol,Q,
                   loose,strict,p_a,mab_p,mab_C,tau,looseDown,looseUp,strictDown,strictUp,p_aDown,p_aUp)


# boxplot PA for varying tau, Q0 ####
tau_type <- c(4,1)
Q_type <- c(10,25)
if (draw_figs){
  par(mfrow=c(2,3))
  colors <- hcl.colors(5,palette="ag_Sunset")
  ylim_mat <- matrix(c(0,1,0,1,0,1,0.4,1,0.4,1,0.4,1),ncol=2,byrow=TRUE)
  k <- 0
  for (distr in distributionType){
    for (sampl in samplingIntensity){
      k <- k+1
      myplot <- boxplot(p_a ~ factor(strategy,levels=strategyType) * factor(tau_sol,levels=tau_type) * factor(Q,levels=Q_type),
                        subset(data, distribution_type == distr & intensity==sampl),
                        col=colors,
                        boxwex=0.8,  notch=TRUE, ylim=ylim_mat[k,],
                        xaxt="n",xlab="")
      
      title(paste(distr,' distribution - intensity ',sampl,sep=""))
      
      my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) paste("tau",x[[2]],"-Q",x[[3]],sep="") )
      my_names <- my_names[seq(1 , length(my_names) , 5)]
      
       axis(1, 
            at = seq(2.5 , 5*length(my_names) , 5), 
            labels = my_names, 
            tick=FALSE , cex=0.3)
       # Add the grey vertical lines
       for(i in seq(0.5 , 21 , 5)){ 
         abline(v=i,lty=1, col="grey")
       }
       for (i in seq(0,1,0.2)){
         abline(h=i,lty=1, col="grey")
       }
    }
  }
}



