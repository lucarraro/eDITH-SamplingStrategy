## eDITH functions
eval_conc <- function(OCN,tau,p,aggregation_level, whichNodes=1:length(p), normalize=FALSE){
  # tau must be in seconds!
  if (aggregation_level=="AG"){
    Nnodes <- OCN$AG$nNodes
    Upstream <- OCN$AG$upstream
    PathLength <- OCN$AG$downstreamPathLength
    AS <- OCN$AG$leng * OCN$AG$width
    velocities <- OCN$AG$pathvelocities
    Q <- OCN$AG$depth * OCN$AG$velocity * OCN$AG$width
  } else if (aggregation_level=="RN"){
    Nnodes <- OCN$RN$nNodes
    Upstream <- OCN$RN$upstream
    PathLength <- as.dgCMatrix.spam(OCN$RN$downstreamPathLength)
    AS <- OCN$RN$leng * OCN$RN$width
    velocities <- as.dgCMatrix.spam(OCN$RN$pathvelocities)
    Q <- OCN$RN$depth * OCN$RN$velocity * OCN$RN$width
  }
  C <- sapply(whichNodes, eDITH, Upstream=Upstream, Q=Q, PathLength=PathLength, 
              velocities=velocities, tau=tau, p=p, AS=AS)
  if (normalize){
    maxC <- sum(AS)/max(Q) 
    C <- C/maxC
  }
  invisible(C)
}

eDITH <- function(i, Upstream, Q, PathLength, velocities, tau, p, AS){
  # tau must be in seconds!
  sources <- Upstream[[i]]
  contrib <- exp(- PathLength[sources,i] / (velocities[sources,i] * tau)) * p[sources] * AS[sources] 
  conc <- sum(contrib)/Q[i]
  invisible(conc)
}

eval_loglik <- function(params, samplingSites, ConcAG, OCN){
  tau <- (params[1] + 1)*3600 # tau can't be lower than 1 h 
  p <- params[-1]
  conc <- eval_conc(OCN, tau, p, "AG", samplingSites, normalize=TRUE)
  loglik <- sum(log(dnorm(conc, ConcAG[samplingSites], 2)))
  if (!isFALSE(loglik < -1e100 )){ # if loglik is NA or -Inf, set random very low value
    loglik <- -(1+runif(1))*1e100
  }
  return(loglik)
}

eval_loglik_taufixed <- function(params, tau_sol, samplingSites, ConcAG, OCN){
  p <- params
  conc <- eval_conc(OCN, tau_sol, p, "AG", samplingSites, normalize=TRUE)
  loglik <- sum(log(dnorm(conc, ConcAG[samplingSites], 2)))
  if (!isFALSE(loglik < -1e100 )){ # if loglik is NA or -Inf, set random very low value
    loglik <- -(1+runif(1))*1e100
  }
  return(loglik)
}