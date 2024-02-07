
prepMahalanobisDistances <- function(dat, coverVariance, treatment){
  # dat should be a matrix with samples in rows and features in columns
  # coverVariance is the percentage of variance that the retained PCs should cover
  # treatment is a vector of sample treatments, corresponding to each row in dat
  
  require(tidyverse)
  
  ############## 1. Calculate the Eigen features from the well-level data 
  PCA = prcomp(dat, center=T, scale=T) # with scaling/centering - 22 components; without scaling - 12; without both - 11
  RotationMatrix = PCA$rotation
  CumProportion = cumsum(PCA$sdev^2)/sum(PCA$sdev^2)
  
  ##############  2. Find the inverse of the covariance matrix
  
  ## Model the data to 
  PC = length(which(CumProportion<coverVariance))+1
  if(PC > dim(dat)[1]){ PC = dim(dat)[1] } 
  RotationMatrix <- RotationMatrix[,1:PC]
  Model = lm(PCA$x[,1:PC] ~ 0 + treatment)
  
  ## 3.b) get covariance matrix
  Cov = estVar(Model)
  
  ## 3.c) inverse
  invCov = solve(Cov)
  
  return(list(RotationMatrix = RotationMatrix, invCov = invCov))
}

computeMahalanobisDistance <- function(dat, RotationMatrix, invCov, treatment, controlID){
  # dat should be a matrix with samples in rows and features in columns
  # RotationMatrix are the PC loadings for the number of PCs that explain coverVariance amount of variability
  # invCov is the inverse covariance matrix computed from the PCs ~ treatment model
  # treatment is a vector with treatment labels for each row in dat
  # controlID is the treatment group that should be considered controls
  
  # get PC scores using loadings and number PCs computed previously
  dat = dat %*% RotationMatrix
  
  # compute the centroid of control samples
  ctrMean = apply(dat[treatment == controlID, ], 2, mean)
  
  # subtract the control centroid from each sample
  Delta = sweep(dat, 2, as.matrix(ctrMean), "-")
  
  # compute the Mahalanobis distance
  D = apply(Delta, 1, function(x) (x %*% invCov %*% x) %>% sqrt() %>% round(., 3))

}
