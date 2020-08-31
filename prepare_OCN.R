rm(list=ls())

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(OCNet)
library(spam)
library(tictoc)

source("eDITH_functions.R")

# prepare OCN ####
# note: recreating OCN_eDNA takes around 1 day
if (!file.exists("OCN_eDNA.rda")){
  set.seed(7)
  OCN <- create_OCN(400, 400, typeInitialState = "T", coolingRate = 0.5, initialNoCoolingPhase = 0.1, cellsize = 50, 
                    displayUpdates = 2)
  OCN <- landscape_OCN(OCN,zMin = 400,slope0=0.005, displayUpdates = 2)
  thr <- find_area_threshold_OCN(OCN, maxReachLength=3000, thrValues = OCN$cellsize^2*(200:2000), displayUpdates = 1)
  OCN <- aggregate_OCN(OCN, thrA=thr$thrValues[which(thr$nNodesAG==201)[1]], maxReachLength = thr$maxReachLength)
  OCN <- paths_OCN(OCN, pathsRN=TRUE, includeDownstreamNode = TRUE)
  OCN <- rivergeometry_OCN(OCN, widthMax = 10)
  
  # calculate velocities
  cat("Calculate AG velocities... \n")
  OCN$AG$pathvelocities <- spam(0,OCN$AG$nNodes,OCN$AG$nNodes)
  set_nodes <- matrix(0,OCN$AG$nNodes^2,2)
  set_values <- numeric(OCN$AG$nNodes^2)
  k <- 1
  for (i in 1:OCN$AG$nNodes){
    for (j in 1:OCN$AG$nNodes){
      path <- OCN$AG$downstreamPath[[i]][[j]]
      if (!is.null(path) && !(i == OCN$AG$outlet && j == OCN$AG$outlet)){
        set_values[k] <- OCN$AG$downstreamPathLength[i,j] / (sum(OCN$AG$leng[path] / OCN$AG$velocity[path]))
        set_nodes[k,] <- c(i,j) 
        k <- k + 1
      } else if (i == OCN$AG$outlet && j == OCN$AG$outlet) {
        set_values[k] <- OCN$AG$velocity[i]
        set_nodes[k,] <- c(i,j) 
        k <- k + 1
      }
    }
  }
  set_values <- set_values[-(k:OCN$AG$nNodes^2)]
  set_nodes <- set_nodes[-(k:OCN$AG$nNodes^2),]
  OCN$AG$pathvelocities[set_nodes] <- set_values
  for (i in 1:OCN$AG$nNodes){
    OCN$AG$pathvelocities[i,i] <- OCN$AG$velocity[i] # patch to correct when length of path is null
  }
  
  cat("Calculate RN velocities... \n")
  OCN$RN$pathvelocities <- spam(0,OCN$RN$nNodes,OCN$RN$nNodes)
  set_nodes <- matrix(0,OCN$RN$nNodes^2,2)
  set_values <- numeric(OCN$RN$nNodes^2)
  k <- 1
  for (i in 1:OCN$RN$nNodes){
    cat(sprintf("\r %.3f",i/OCN$RN$nNodes))
    for (j in 1:OCN$RN$nNodes){
      path <- OCN$RN$downstreamPath[[i]][[j]]
      if (!is.null(path) && !(i == OCN$RN$outlet && j == OCN$RN$outlet)){
        set_values[k] <- OCN$RN$downstreamPathLength[i,j] / (sum(OCN$RN$leng[path] / OCN$RN$velocity[path]))
        set_nodes[k,] <- c(i,j) 
        k <- k + 1
      } else if (i == OCN$RN$outlet && j == OCN$RN$outlet) {
        set_values[k] <- OCN$RN$velocity[i]
        set_nodes[k,] <- c(i,j) 
        k <- k + 1
      }
    }
  }
  set_values <- set_values[-(k:OCN$RN$nNodes^2)]
  set_nodes <- set_nodes[-(k:OCN$RN$nNodes^2),]
  OCN$RN$pathvelocities[set_nodes] <- set_values
  for (i in 1:OCN$RN$nNodes){
    OCN$RN$pathvelocities[i,i] <- OCN$RN$velocity[i] # patch to correct when length of path is null
  }
  
  save(OCN,file="OCN_eDNA.rda")
} else {
  load(file="OCN_eDNA.rda")
}


# define colormap for concentrations
jet <- colorRampPalette(c("#0000ff","#0080ff","#00ffff","#80ff80","#ffff00","#ff8000","#ff0000"))

# pointsource RN
theme_pointRN <- numeric(OCN$RN$nNodes)
site <- which(OCN$RN$downstreamPathLength[,OCN$RN$outlet]==max(OCN$RN$downstreamPathLength[,OCN$RN$outlet]))
theme_pointRN[site] <- sum(OCN$RN$leng * OCN$RN$width)/(OCN$RN$leng[site] * OCN$RN$width[site])


# eval conc patterns
Conc_1h_unif <- eval_conc(OCN,3600*1,1+numeric(OCN$RN$nNodes),"RN",normalize=TRUE)
Conc_4h_unif <- eval_conc(OCN,3600*4,1+numeric(OCN$RN$nNodes),"RN",normalize=TRUE)
Conc_Inf_unif <- eval_conc(OCN,Inf,1+numeric(OCN$RN$nNodes),"RN",normalize=TRUE)

Conc_1h_pointRN <- eval_conc(OCN,3600*1,theme_pointRN,"RN",normalize=TRUE)
Conc_4h_pointRN <- eval_conc(OCN,3600*4,theme_pointRN,"RN",normalize=TRUE)
Conc_Inf_pointRN <- eval_conc(OCN,Inf,theme_pointRN,"RN",normalize=TRUE)


# semilog plot Conc vs dist to outlet ####
colors <- hcl.colors(7,palette="RdYlBu")
site <- which(OCN$RN$downstreamPathLength[,OCN$RN$outlet]==max(OCN$RN$downstreamPathLength[,OCN$RN$outlet]))
DtoO <- OCN$RN$downstreamPathLength[OCN$RN$downstreamPath[[site]][[OCN$RN$outlet]],OCN$RN$outlet]
par(mai=c(1,1,1,1), mfrow=c(1,1))
plot(DtoO/1000,log10(Conc_Inf_unif[OCN$RN$downstreamPath[[site]][[OCN$RN$outlet]]]),
     pch=20,xlim=c(40,0),ylim=c(-2,3),axes=FALSE,col=colors[1],type="l",
     xlab="Distance to Outlet [km]",ylab="Relative eDNA concentration [-]")
axis(1,at=c(40,30,20,10,0),pos=-2)
axis(2,at=c(-2,-1,0,1,2,3),pos=40)
points(DtoO/1000,log10(Conc_4h_unif[OCN$RN$downstreamPath[[site]][[OCN$RN$outlet]]]),
       pch=20,col=colors[2],type="l")
points(DtoO/1000,log10(Conc_1h_unif[OCN$RN$downstreamPath[[site]][[OCN$RN$outlet]]]),
       pch=20,col=colors[3],type="l")

points(DtoO/1000,log10(Conc_Inf_pointRN[OCN$RN$downstreamPath[[site]][[OCN$RN$outlet]]]),
       pch=20,col=colors[7],type="l")
points(DtoO/1000,log10(Conc_4h_pointRN[OCN$RN$downstreamPath[[site]][[OCN$RN$outlet]]]),
       pch=20,col=colors[6],type="l")
points(DtoO/1000,log10(Conc_1h_pointRN[OCN$RN$downstreamPath[[site]][[OCN$RN$outlet]]]),
       pch=20,col=colors[5],type="l")



## concentration maps ####
par(mfrow=c(3,2), mai=c(0,0,0,0))
draw_thematic_OCN(log10(Conc_Inf_unif),OCN,colLevels=c(-1,1,1000),colPalette=jet)
draw_thematic_OCN(log10(Conc_Inf_pointRN),OCN,colLevels=c(-1,1,1000),colPalette=jet)
draw_thematic_OCN(log10(Conc_4h_unif),OCN,colLevels=c(-1,1,1000),colPalette=jet)
draw_thematic_OCN(log10(Conc_4h_pointRN),OCN,colLevels=c(-1,1,1000),colPalette=jet)
draw_thematic_OCN(log10(Conc_1h_unif),OCN,colLevels=c(-1,1,1000),colPalette=jet)
draw_thematic_OCN(log10(Conc_1h_pointRN),OCN,colLevels=c(-1,1,1000),colPalette=jet)

# initial figure ####
par(mfrow=c(1,1))
Nodes <- setdiff(1:OCN$AG$nNodes,OCN$AG$outlet)
NodesDownstream <- intersect(which(OCN$AG$A > median(OCN$AG$A[Nodes])), Nodes)
NodesUpstream <- intersect(which(OCN$AG$A <= median(OCN$AG$A[Nodes])), Nodes)
theme <- numeric(OCN$AG$nNodes)
theme[NodesDownstream] <- 1
draw_thematic_OCN(theme,OCN,drawNodes=TRUE,nodeType="downstream",cex=1)
