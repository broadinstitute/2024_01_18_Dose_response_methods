# Compute tPODs for all BMD lists
# Jessica Ewald
# February 2, 2024


# The purpose of this script is to compute these types of tPODs for each list
# of BMDs: 10th percentile, first mode, LCRD, pathway, global mahalanobis, and 
# PLS (new). 

#### Setup ####

library(dplyr)
library(arrow)
source("../fastbmdR/fastbmdR_main.R")
source("../fastbmdR/fastbmdR_utils.R")

# relative path to data submodule
data.path <- "../2024_01_18_Dose_response_methods-data/data/"

# input filepaths
bmd.paths <- c(paste0(data.path, "4_bmd_results/httr_bmd_aov.parquet"),
               paste0(data.path, "4_bmd_results/httr_bmd_wtt.parquet"),
               paste0(data.path, "4_bmd_results/httr_bmd_s1500.parquet"),
               paste0(data.path, "4_bmd_results/httr_bmd_nomic.parquet"),
               paste0(data.path, "4_bmd_results/htpp_bmd_all.parquet"),
               paste0(data.path, "4_bmd_results/htpp_bmd_aov.parquet"))

names(bmd.paths) <- c("httr_aov", "httr_wtt", "httr_s1500", "httr_nomic", "htpp_all", "htpp_aov")


#### Compute tPODs ####

tpod.res <- data.frame()
for(i in c(1:length(bmd.paths))){
  
  # compute 10th centile
  
  # compute mode
  
  # compute lcrd
  
  # compute pathway (go bp - keep other gene set analysis separate)
  if(grepl("httr", names(bmd.paths[i]))){
    
  }
  
  # compute global mahalanobis
  
  # compute pls
  
}


#### Compute differences in tPOD distributions across conditions ####

# Metrics: 
# KS p-value for tPOD dists between conditions, for each of the 11 chemicals (paired?)
# Paired difference of tPOD dist mean (1 test for all 11 chemicals) between conditions
# Paired difference of tPOD dist variance (1 test for all 11 chemicals) between conditions

# Across tPOD types for the same dataset

# Across datasets for the same tpod type

