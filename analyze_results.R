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

simulationType <- c("tf","tu","me") # 
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


load("results_all.rda")

simulation <- vector("character",nTotSim)
strategy <- distribution_type <- distribution_no <- simulation
intensity <- vector("numeric",nTotSim)
loose <- strict <- strict2 <- strict3 <- p_a <- p_a2 <- p_a3 <- mab_p <- mab_C <- tau <- intensity
looseDown <- looseUp <- strictDown <- strictUp  <- p_aDown <- p_aUp <- intensity

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
              
              strict2[k] <- sum( (p_sol <= 0.5 & p <= 0.5) | (p >= 0.5*p_sol & p <= 2*p_sol), na.rm=T) / 200
              strict3[k] <- sum( (p_sol <= 0.125 & p <= 0.125) | (p >= 9/10*p_sol & p <= 10/9*p_sol), na.rm=T) / 200
              strict[k] <- sum( (p_sol <= 0.25 & p <= 0.25) | (p >= 0.75*p_sol & p <= 4/3*p_sol), na.rm=T) / 200
              p_a[k] <- sum( (p_sol <= 0.25 & p <= 0.25) | (p_sol > 0.2 & p > 0.2), na.rm=T) / 200
              p_a2[k] <- sum( (p_sol <= 0.5 & p <= 0.5) | (p_sol > 0.45 & p > 0.45), na.rm=T) / 200
              p_a3[k] <- sum( (p_sol <= 0.125 & p <= 0.125) | (p_sol > 0.075 & p > 0.075), na.rm=T) / 200
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
data <- data.frame(simulation,strategy,intensity,distribution_type,distribution_no,strict2,strict3,
                   strict,p_a,p_a2,p_a3,mab_p,mab_C,tau,looseDown,looseUp,strictDown,strictUp,p_aDown,p_aUp)


# boxplot PA ####
if (draw_figs){
  par(mfrow=c(2,3))
  colors <- hcl.colors(5,palette="ag_Sunset")
  ylim_mat <- matrix(c(0,1,0,1,0,1,0.4,1,0.4,1,0.4,1),ncol=2,byrow=TRUE)
  k <- 0
  for (distr in distributionType){
    for (sampl in samplingIntensity){
      k <- k+1
      myplot <- boxplot(p_a ~ factor(strategy,levels=strategyType) * factor(simulation,levels=simulationType),
                        subset(data, distribution_type == distr & intensity==sampl),
                        col=colors,
                        boxwex=0.8,  notch=TRUE, ylim=ylim_mat[k,],
                        xaxt="n",xlab="")
      
      title(paste(distr,' distribution - intensity ',sampl,sep=""))
      
      my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
      my_names <- my_names[seq(1 , length(my_names) , 5)]
      
      axis(1, 
           at = seq(2.5 , 5*length(my_names) , 5), 
           labels = my_names, 
           tick=FALSE , cex=0.3)
      # Add the grey vertical lines
      for(i in seq(0.5 , 20 , 5)){ 
        abline(v=i,lty=1, col="grey")
      }
      for (i in seq(0,1,0.2)){
        abline(h=i,lty=1, col="grey")
      }
    }
  }
}

# boxplot D ####
if(draw_figs){
  par(mfrow=c(2,3))
  colors <- hcl.colors(5,palette="ag_Sunset")
  ylim_mat <- matrix(c(0,1,0,1,0,1,0,0.8,0,0.8,0,0.8),ncol=2,byrow=TRUE)
  k <- 0
  for (distr in distributionType){
    for (sampl in samplingIntensity){
      k <- k+1
      myplot <- boxplot(strict ~ factor(strategy,levels=strategyType) * factor(simulation,levels=simulationType),
                        subset(data, distribution_type == distr & intensity==sampl),
                        col=colors,
                        boxwex=0.8,  notch=TRUE, ylim=ylim_mat[k,],
                        xaxt="n",xlab="")
      
      title(paste(distr,' distribution - intensity ',sampl,sep=""))
      
      my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
      my_names <- my_names[seq(1 , length(my_names) , 5)]
      
      axis(1, 
           at = seq(2.5 , 5*length(my_names) , 5), 
           labels = my_names, 
           tick=FALSE , cex=0.3)
      # Add the grey vertical lines
      for(i in seq(0.5 , 20 , 5)){ 
        abline(v=i,lty=1, col="grey")
      }
      for (i in seq(0,1,0.2)){
        abline(h=i,lty=1, col="grey")
      }
    }
  }
}

# boxplot tau ####
if(draw_figs){
  par(mfrow=c(2,3))
  colors <- hcl.colors(5,palette="ag_Sunset")
  #ylim_mat <- matrix(c(0,0.8,0,0.4,0.2,1,0,0.6,0.4,1,0,0.8),ncol=2,byrow=TRUE)
  k <- 0
  
  for (distr in distributionType){
    for (sampl in samplingIntensity){
      k <- k+1
      myplot <- boxplot(tau ~ factor(strategy,levels=strategyType) * factor(simulation,levels=simulationType),
                        subset(data, distribution_type == distr & intensity==sampl),
                        col=colors,
                        boxwex=0.8,  #notch=TRUE, #ylim=c(0,1),#ylim_mat[k,],
                        xaxt="n",xlab="",log="y")
      
      title(paste(distr,' distribution - intensity ',sampl,sep=""))
      
      c
      # Add the grey vertical lines
      for(i in seq(0.5 , 20 , 5)){
        abline(v=i,lty=1, col="grey")
      }
      # for (i in seq(0,1,0.2)){
      #   abline(h=i,lty=1, col="grey")
      # }
      # Add a legend
      if (k==6){
        legend("bottomright", legend = strategyType,
               col=colors,
               pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = c(0.1, 0.1))
      }
    }
  }
}

# boxplot mab_p ####
if(draw_figs){
  par(mfrow=c(2,3))
  colors <- hcl.colors(5,palette="ag_Sunset")
  ylim_mat <- matrix(c(-0.5,1.25,-0.5,1.25,-0.5,1.25,-0.5,1.25,-0.5,1.25,-0.5,1.25),ncol=2,byrow=TRUE)
  k <- 0
  for (distr in distributionType){
    for (sampl in samplingIntensity){
      k <- k+1
      myplot <- boxplot(log10(mab_p) ~ factor(strategy,levels=strategyType) * factor(simulation,levels=simulationType),
                        subset(data, distribution_type == distr & intensity==sampl),
                        col=colors,
                        boxwex=0.8,  ylim=ylim_mat[k,],
                        xaxt="n",xlab="")
      
      title(paste(distr,' distribution - intensity ',sampl,sep=""))
      
      my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
      my_names <- my_names[seq(1 , length(my_names) , 5)]
      
      axis(1, 
           at = seq(2.5 , 5*length(my_names) , 5), 
           labels = my_names, 
           tick=FALSE , cex=0.3)
      # Add the grey vertical lines
      for(i in seq(0.5 , 20 , 5)){ 
        abline(v=i,lty=1, col="grey")
      }
      for (i in seq(-0.5,1,0.5)){
        abline(h=i,lty=1, col="grey")
      }
    }
  }
}

# boxplot PA2 ####
if (draw_figs){
  par(mfrow=c(2,3))
  colors <- hcl.colors(5,palette="ag_Sunset")
  ylim_mat <- matrix(c(0,1,0,1,0,1,0.4,1,0.4,1,0.4,1),ncol=2,byrow=TRUE)
  k <- 0
  for (distr in distributionType){
    for (sampl in samplingIntensity){
      k <- k+1
      myplot <- boxplot(p_a2 ~ factor(strategy,levels=strategyType) * factor(simulation,levels=simulationType),
                        subset(data, distribution_type == distr & intensity==sampl),
                        col=colors,
                        boxwex=0.8,  notch=TRUE, ylim=ylim_mat[k,],
                        xaxt="n",xlab="")
      
      title(paste(distr,' distribution - intensity ',sampl,sep=""))
      
      my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
      my_names <- my_names[seq(1 , length(my_names) , 5)]
      
      axis(1, 
           at = seq(2.5 , 5*length(my_names) , 5), 
           labels = my_names, 
           tick=FALSE , cex=0.3)
      # Add the grey vertical lines
      for(i in seq(0.5 , 20 , 5)){ 
        abline(v=i,lty=1, col="grey")
      }
      for (i in seq(0,1,0.2)){
        abline(h=i,lty=1, col="grey")
      }
    }
  }
}

# boxplot PA3 ####
if (draw_figs){
  par(mfrow=c(2,3))
  colors <- hcl.colors(5,palette="ag_Sunset")
  ylim_mat <- matrix(c(0,1,0,1,0,1,0.4,1,0.4,1,0.4,1),ncol=2,byrow=TRUE)
  k <- 0
  for (distr in distributionType){
    for (sampl in samplingIntensity){
      k <- k+1
      myplot <- boxplot(p_a3 ~ factor(strategy,levels=strategyType) * factor(simulation,levels=simulationType),
                        subset(data, distribution_type == distr & intensity==sampl),
                        col=colors,
                        boxwex=0.8,  notch=TRUE, ylim=ylim_mat[k,],
                        xaxt="n",xlab="")
      
      title(paste(distr,' distribution - intensity ',sampl,sep=""))
      
      my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
      my_names <- my_names[seq(1 , length(my_names) , 5)]
      
      axis(1, 
           at = seq(2.5 , 5*length(my_names) , 5), 
           labels = my_names, 
           tick=FALSE , cex=0.3)
      # Add the grey vertical lines
      for(i in seq(0.5 , 20 , 5)){ 
        abline(v=i,lty=1, col="grey")
      }
      for (i in seq(0,1,0.2)){
        abline(h=i,lty=1, col="grey")
      }
    }
  }
}

# boxplot D2 ####
if(draw_figs){
  par(mfrow=c(2,3))
  colors <- hcl.colors(5,palette="ag_Sunset")
  ylim_mat <- matrix(c(0,1,0,1,0,1,0,0.8,0,0.8,0,0.8),ncol=2,byrow=TRUE)
  k <- 0
  for (distr in distributionType){
    for (sampl in samplingIntensity){
      k <- k+1
      myplot <- boxplot(strict2 ~ factor(strategy,levels=strategyType) * factor(simulation,levels=simulationType),
                        subset(data, distribution_type == distr & intensity==sampl),
                        col=colors,
                        boxwex=0.8,  notch=TRUE, ylim=ylim_mat[k,],
                        xaxt="n",xlab="")
      
      title(paste(distr,' distribution - intensity ',sampl,sep=""))
      
      my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
      my_names <- my_names[seq(1 , length(my_names) , 5)]
      
      axis(1, 
           at = seq(2.5 , 5*length(my_names) , 5), 
           labels = my_names, 
           tick=FALSE , cex=0.3)
      # Add the grey vertical lines
      for(i in seq(0.5 , 20 , 5)){ 
        abline(v=i,lty=1, col="grey")
      }
      for (i in seq(0,1,0.2)){
        abline(h=i,lty=1, col="grey")
      }
    }
  }
}

# boxplot D3 ####
if(draw_figs){
  par(mfrow=c(2,3))
  colors <- hcl.colors(5,palette="ag_Sunset")
  ylim_mat <- matrix(c(0,1,0,1,0,1,0,0.8,0,0.8,0,0.8),ncol=2,byrow=TRUE)
  k <- 0
  for (distr in distributionType){
    for (sampl in samplingIntensity){
      k <- k+1
      myplot <- boxplot(strict3 ~ factor(strategy,levels=strategyType) * factor(simulation,levels=simulationType),
                        subset(data, distribution_type == distr & intensity==sampl),
                        col=colors,
                        boxwex=0.8,  notch=TRUE, ylim=ylim_mat[k,],
                        xaxt="n",xlab="")
      
      title(paste(distr,' distribution - intensity ',sampl,sep=""))
      
      my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
      my_names <- my_names[seq(1 , length(my_names) , 5)]
      
      axis(1, 
           at = seq(2.5 , 5*length(my_names) , 5), 
           labels = my_names, 
           tick=FALSE , cex=0.3)
      # Add the grey vertical lines
      for(i in seq(0.5 , 20 , 5)){ 
        abline(v=i,lty=1, col="grey")
      }
      for (i in seq(0,1,0.2)){
        abline(h=i,lty=1, col="grey")
      }
    }
  }
}

# plot distribution types ####
if(draw_figs){
  par(mfrow=c(2,5))
  draw_thematic_OCN(results_all$tf$u0$s10$R1$`1`$p_sol, OCN, colLevels=c(0,10,1000)); title("R1")
  draw_thematic_OCN(results_all$tf$u0$s10$R2$`1`$p_sol, OCN, colLevels=c(0,10,1000)); title("R2")
  draw_thematic_OCN(results_all$tf$u0$s10$R3$`1`$p_sol, OCN, colLevels=c(0,10,1000)); title("R3")
  draw_thematic_OCN(results_all$tf$u0$s10$R4$`1`$p_sol, OCN, colLevels=c(0,10,1000)); title("R4")
  draw_thematic_OCN(results_all$tf$u0$s10$R5$`1`$p_sol, OCN, colLevels=c(0,10,1000)); title("R5")
  draw_thematic_OCN(results_all$tf$u0$s10$A1$`1`$p_sol, OCN, colLevels=c(0,10,1000)); title("A1")
  draw_thematic_OCN(results_all$tf$u0$s10$A2$`1`$p_sol, OCN, colLevels=c(0,10,1000)); title("A2")
  draw_thematic_OCN(results_all$tf$u0$s10$A3$`1`$p_sol, OCN, colLevels=c(0,10,1000)); title("A3")
  draw_thematic_OCN(results_all$tf$u0$s10$A4$`1`$p_sol, OCN, colLevels=c(0,10,1000)); title("A4")
  draw_thematic_OCN(results_all$tf$u0$s10$A5$`1`$p_sol, OCN, colLevels=c(0,10,1000)); title("A5")
}

# plot C_sim maps ####
if(draw_figs){
  jet <- colorRampPalette(c("#0000ff","#0080ff","#00ffff","#80ff80","#ffff00","#ff8000","#ff0000"))
  par(mfrow=c(2,5))
  draw_thematic_OCN(log10(results_all$tf$u0$s10$R1$`1`$C_sol), OCN, colPalette=jet, colLevels=c(-1,1,1000)); title("R1")
  draw_thematic_OCN(log10(results_all$tf$u0$s10$R2$`1`$C_sol), OCN, colPalette=jet, colLevels=c(-1,1,1000)); title("R2")
  draw_thematic_OCN(log10(results_all$tf$u0$s10$R3$`1`$C_sol), OCN, colPalette=jet, colLevels=c(-1,1,1000)); title("R3")
  draw_thematic_OCN(log10(results_all$tf$u0$s10$R4$`1`$C_sol), OCN, colPalette=jet, colLevels=c(-1,1,1000)); title("R4")
  draw_thematic_OCN(log10(results_all$tf$u0$s10$R5$`1`$C_sol), OCN, colPalette=jet, colLevels=c(-1,1,1000)); title("R5")
  draw_thematic_OCN(log10(results_all$tf$u0$s10$A1$`1`$C_sol), OCN, colPalette=jet, colLevels=c(-1,1,1000)); title("A1")
  draw_thematic_OCN(log10(results_all$tf$u0$s10$A2$`1`$C_sol), OCN, colPalette=jet, colLevels=c(-1,1,1000)); title("A2")
  draw_thematic_OCN(log10(results_all$tf$u0$s10$A3$`1`$C_sol), OCN, colPalette=jet, colLevels=c(-1,1,1000)); title("A3")
  draw_thematic_OCN(log10(results_all$tf$u0$s10$A4$`1`$C_sol), OCN, colPalette=jet, colLevels=c(-1,1,1000)); title("A4")
  draw_thematic_OCN(log10(results_all$tf$u0$s10$A5$`1`$C_sol), OCN, colPalette=jet, colLevels=c(-1,1,1000)); title("A5")
}

# boxplot on nestedness ####

# pick sampling sites
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


# calculate distance matrix
dist_conn <- OCN$AG$downstreamPathLength
dist_unconn <- OCN$AG$downstreamLengthUnconnected

dist_conn <- dist_conn + t(dist_conn) - diag(diag(dist_conn)) # make it symmetric, include all connected sites either downstream or upstream
dist_unconn <- dist_unconn[dist_conn==0] # remove upstream connected sites
dist_unconn <- dist_unconn + t(dist_unconn) # make it symmetric

distanz <- dist_conn + dist_unconn


# calculate nestedness and mean distance for each set of sampling sites
nestedness <- mean_distanz <- numeric(10*5*3)
sampling <- strategy <- vector("character",10*5*3)
k<-0
for (sampl in samplingIntensity){
  for (strat in strategyType){
    for (indSim in 1:10){
      k <- k+1
      tmp <- samplingSitesList[[sampl]][[strat]][indSim,]
      for (i in tmp){
        path <-  OCN$AG$downstreamPath[[i]][[outlet]]
        nestedness[k] <- nestedness[k] + length(setdiff(intersect(path,tmp),i))
      }
      nestedness[k] <- nestedness[k]/(sum(1:(length(tmp)-1)))
      sampling[k] <- sampl
      strategy[k] <- strat
      
      # mean distance between sites
      all_pairs <- combn(tmp, 2) # all possible pairs of sites (counted only once!)
      mean_dist <- 0
      for (i in 1:(length(all_pairs)/2)){
        mean_dist <- mean_dist + distanz[all_pairs[1,i],all_pairs[2,i]]
      }
      mean_dist <- mean_dist/(length(all_pairs)/2)
      mean_distanz[k] <- mean_dist
    }
  }
}
data_conn <- data.frame(nestedness, mean_distanz, sampling, strategy)

# boxplot nestedness
if(draw_figs){
  myplot <- boxplot(nestedness ~ factor(strategy,levels=strategyType)*factor(sampling,levels=samplingIntensity),
                    col=colors, boxwex=0.8,xaxt="n",xlab="Sampling intensity",ylim=c(0,0.4))
  my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
  my_names <- my_names[seq(1 , length(my_names) , 5)]
  
  axis(1,
       at = seq(2.5 , 5*length(my_names) , 5),
       labels = my_names,
       tick=FALSE , cex=0.3)
  for(i in seq(0.5 , 20 , 5)){ 
    abline(v=i,lty=1, col="grey")
  }
}


# boxplot mean_distanz
if(draw_figs){
  myplot <- boxplot(mean_distanz ~ factor(strategy,levels=strategyType)*factor(sampling,levels=samplingIntensity),
                    col=colors, boxwex=0.8,xaxt="n",xlab="Sampling intensity",ylim=c(20000,35000))
  my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
  my_names <- my_names[seq(1 , length(my_names) , 5)]
  
  axis(1,
       at = seq(2.5 , 5*length(my_names) , 5),
       labels = my_names,
       tick=FALSE , cex=0.3)
  for(i in seq(0.5 , 20 , 5)){ 
    abline(v=i,lty=1, col="grey")
  }
}


# plot examples of S maps ####
if(draw_figs){
  #pdf("example_S.pdf",width=16,height=8,paper="special")
  par(mfrow=c(2,4))
  draw_thematic_OCN(results_all$me$u0$s10$R1$`1`$p_sol, OCN, colLevels=c(0,10,1000)); title("R1")
  draw_thematic_OCN(results_all$me$u20$s10$R1$`1`$par[-1], OCN, colLevels=c(0,10,1000),addLegend = FALSE); 
  title(paste("S10 / U20 / PA = ",results_all$me$u20$s10$R1$`1`$p_a," / D = ",results_all$me$u20$s10$R1$`1`$strict,sep=""))
  points(OCN$AG$XReach[results_all$me$u20$s10$R1$`1`$samplingSites],
         OCN$AG$YReach[results_all$me$u20$s10$R1$`1`$samplingSites],col="blue")
  draw_thematic_OCN(results_all$me$u20$s25$R1$`1`$par[-1], OCN, colLevels=c(0,10,1000),addLegend = FALSE); 
  title(paste("S25 / U20 / PA = ",results_all$me$u20$s25$R1$`1`$p_a," / D = ",results_all$me$u20$s25$R1$`1`$strict,sep=""))
  points(OCN$AG$XReach[results_all$me$u20$s25$R1$`1`$samplingSites],
         OCN$AG$YReach[results_all$me$u20$s25$R1$`1`$samplingSites],col="blue")
  draw_thematic_OCN(results_all$me$u20$s50$R1$`1`$par[-1], OCN, colLevels=c(0,10,1000),addLegend = FALSE); 
  title(paste("S50 / U20 / PA = ",results_all$me$u20$s50$R1$`1`$p_a," / D = ",results_all$me$u20$s50$R1$`1`$strict,sep=""))
  points(OCN$AG$XReach[results_all$me$u20$s50$R1$`1`$samplingSites],
         OCN$AG$YReach[results_all$me$u20$s50$R1$`1`$samplingSites],col="blue")
  frame()
  draw_thematic_OCN(results_all$me$u80$s10$R1$`1`$par[-1], OCN, colLevels=c(0,10,1000),addLegend = FALSE); 
  title(paste("S10 / U80 / PA = ",results_all$me$u80$s10$R1$`1`$p_a," / D = ",results_all$me$u80$s10$R1$`1`$strict,sep=""))
  points(OCN$AG$XReach[results_all$me$u80$s10$R1$`1`$samplingSites],
         OCN$AG$YReach[results_all$me$u80$s10$R1$`1`$samplingSites],col="blue")
  draw_thematic_OCN(results_all$me$u80$s25$R1$`1`$par[-1], OCN, colLevels=c(0,10,1000),addLegend = FALSE); 
  title(paste("S25 / U80 / PA = ",results_all$me$u80$s25$R1$`1`$p_a," / D = ",results_all$me$u80$s25$R1$`1`$strict,sep=""))
  points(OCN$AG$XReach[results_all$me$u80$s25$R1$`1`$samplingSites],
         OCN$AG$YReach[results_all$me$u80$s25$R1$`1`$samplingSites],col="blue")
  draw_thematic_OCN(results_all$me$u80$s50$R1$`1`$par[-1], OCN, colLevels=c(0,10,1000),addLegend = FALSE); 
  title(paste("S50 / U80 / PA = ",results_all$me$u80$s50$R1$`1`$p_a," / D = ",results_all$me$u80$s50$R1$`1`$strict,sep=""))
  points(OCN$AG$XReach[results_all$me$u80$s50$R1$`1`$samplingSites],
         OCN$AG$YReach[results_all$me$u80$s50$R1$`1`$samplingSites],col="blue")
  #dev.off()
}

# plot examples of E maps ####
if(draw_figs){
  #pdf("example_E.pdf",width=16,height=8,paper="special")
  par(mfrow=c(2,4))
  draw_thematic_OCN(results_all$me$u0$s10$A1$`1`$p_sol, OCN, colLevels=c(0,10,1000)); title("A1")
  draw_thematic_OCN(results_all$me$u20$s10$A1$`1`$par[-1], OCN, colLevels=c(0,10,1000),addLegend = FALSE); 
  title(paste("S10 / U20 / PA = ",results_all$me$u20$s10$A1$`1`$p_a," / D = ",results_all$me$u20$s10$A1$`1`$strict,sep=""))
  points(OCN$AG$XReach[results_all$me$u20$s10$A1$`1`$samplingSites],
         OCN$AG$YReach[results_all$me$u20$s10$A1$`1`$samplingSites],col="blue")
  draw_thematic_OCN(results_all$me$u20$s25$A1$`1`$par[-1], OCN, colLevels=c(0,10,1000),addLegend = FALSE); 
  title(paste("S25 / U20 / PA = ",results_all$me$u20$s25$A1$`1`$p_a," / D = ",results_all$me$u20$s25$A1$`1`$strict,sep=""))
  points(OCN$AG$XReach[results_all$me$u20$s25$A1$`1`$samplingSites],
         OCN$AG$YReach[results_all$me$u20$s25$A1$`1`$samplingSites],col="blue")
  draw_thematic_OCN(results_all$me$u20$s50$A1$`1`$par[-1], OCN, colLevels=c(0,10,1000),addLegend = FALSE); 
  title(paste("S50 / U20 / PA = ",results_all$me$u20$s50$A1$`1`$p_a," / D = ",results_all$me$u20$s50$A1$`1`$strict,sep=""))
  points(OCN$AG$XReach[results_all$me$u20$s50$A1$`1`$samplingSites],
         OCN$AG$YReach[results_all$me$u20$s50$A1$`1`$samplingSites],col="blue")
  frame()
  draw_thematic_OCN(results_all$me$u80$s10$A1$`1`$par[-1], OCN, colLevels=c(0,10,1000),addLegend = FALSE); 
  title(paste("S10 / U80 / PA = ",results_all$me$u80$s10$A1$`1`$p_a," / D = ",results_all$me$u80$s10$A1$`1`$strict,sep=""))
  points(OCN$AG$XReach[results_all$me$u80$s10$A1$`1`$samplingSites],
         OCN$AG$YReach[results_all$me$u80$s10$A1$`1`$samplingSites],col="blue")
  draw_thematic_OCN(results_all$me$u80$s25$A1$`1`$par[-1], OCN, colLevels=c(0,10,1000),addLegend = FALSE); 
  title(paste("S25 / U80 / PA = ",results_all$me$u80$s25$A1$`1`$p_a," / D = ",results_all$me$u80$s25$A1$`1`$strict,sep=""))
  points(OCN$AG$XReach[results_all$me$u80$s25$A1$`1`$samplingSites],
         OCN$AG$YReach[results_all$me$u80$s25$A1$`1`$samplingSites],col="blue")
  draw_thematic_OCN(results_all$me$u80$s50$A1$`1`$par[-1], OCN, colLevels=c(0,10,1000),addLegend = FALSE); 
  title(paste("S50 / U80 / PA = ",results_all$me$u80$s50$A1$`1`$p_a," / D = ",results_all$me$u80$s50$A1$`1`$strict,sep=""))
  points(OCN$AG$XReach[results_all$me$u80$s50$A1$`1`$samplingSites],
         OCN$AG$YReach[results_all$me$u80$s50$A1$`1`$samplingSites],col="blue")
  #dev.off()
}


