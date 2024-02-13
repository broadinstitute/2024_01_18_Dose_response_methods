# Compute tPODs for all BMD lists
# Jessica Ewald
# February 2, 2024


# The purpose of this script is to compute these types of tPODs for each list
# of BMDs: 10th percentile, first mode, LCRD, pathway, global mahalanobis, and 
# PLS (new). 

#### Setup ####

library(dplyr)
library(arrow)
library(edgeR)
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
               paste0(data.path, "4_bmd_results/htpp_bmd_aov.parquet"),
               paste0(data.path, "4_bmd_results/htpp_bmd_wtt.parquet"))

mat.paths <- c(list.files(paste0(data.path, "2_sampled_data/httr/2_reps"), full.names = T),
               list.files(paste0(data.path, "2_sampled_data/httr/3_reps"), full.names = T),
               list.files(paste0(data.path, "2_sampled_data/htpp/2_reps"), full.names = T),
               list.files(paste0(data.path, "2_sampled_data/htpp/3_reps"), full.names = T))

names(bmd.paths) <- c("httr_aov", "httr_wtt", "httr_s1500", "httr_nomic", "htpp_all", "htpp_aov", "htpp_wtt")
feature.set <- c("ANOVA", "WTT", "S1500", "Nomic", "All", "ANOVA", "WTT")

# get gene set info
probe.map <- read.csv("/Users/jessicaewald/NetbeansProjects/2024_01_18_Dose_response_methods/2024_01_18_Dose_response_methods-data/data/probemaps/100739_homo_sapiens_wt_probemap.csv")
universe <- probe.map$PROBE_NAME %>% unique()
gs.list <- readRDS("/Users/jessicaewald/NetbeansProjects/2024_01_18_Dose_response_methods/2024_01_18_Dose_response_methods-data/data/gs_libraries/go_bp.rds")[["sets"]]

# get global MD pre-reqs
# full httr matrix
httr <- read_parquet("../2024_01_18_Dose_response_methods-data/data/1_orig_data/httr_counts.parquet") %>% as.data.frame()
httr.meta <- read_parquet("../2024_01_18_Dose_response_methods-data/data/1_orig_data/httr_key.parquet") %>% as.data.frame()
rownames(httr) <- httr$probe_id
httr <- httr[,-1]
httr <- as.matrix(httr)
httr <- httr[,colnames(httr) %in% httr.meta$sample_id[httr.meta$qc_flag == "OK"]]
httr.meta <- httr.meta[httr.meta$qc_flag == "OK", ]

# filter based on abundance
keep.feat <- apply(httr, 1, mean) > 5
httr <- httr[keep.feat, ]

# convert to log-cpm
nf <- calcNormFactors(httr, method="RLE")
httr <- cpm(httr, lib.size=colSums(httr)*nf, log = TRUE, prior.count = 1)
httr <- t(httr)

# Check out htpp data
htpp <- read_parquet("../2024_01_18_Dose_response_methods-data/data/1_orig_data/htpp_well_norm.parquet") %>% as.data.frame()
htpp <- htpp[htpp$plate_id == "TC00000483", ]
htpp.meta <- htpp[,c(1:11)]
htpp <- htpp[,-c(1:11)]
rownames(htpp) <- htpp.meta$sample_id
htpp <- htpp[,!(colnames(htpp) %in% c("f_1073", "f_1074"))]

# compute MD stats
httr.gmd <- prepMahalanobisDistances(httr, 0.95, paste0(httr.meta$chem_name, "_", httr.meta$dose_level))
htpp.gmd <- prepMahalanobisDistances(htpp, 0.95, paste0(htpp.meta$chem_id, "_", htpp.meta$dose_level))

#### Compute tPODs ####

tpod.res <- data.frame()
httr.pls.loadings <- list()
htpp.pls.loadings <- list()
for(i in c(1:length(bmd.paths))){
  bmds <- read_parquet(bmd.paths[i])
  print(bmd.paths[i])
  
  analysisIDs <- unique(bmds$analysisID)
  
  info <- lapply(analysisIDs, function(x) {unlist(strsplit(x, "_"))})
  info <- unlist(info)
  info <- matrix(info, ncol = 4, byrow = TRUE) %>% as.data.frame()
  colnames(info) <- c("Compound", "Replicates", "Iteration", "Platform")
  info$Features <- feature.set[i]
  info$AnalysisID <- analysisIDs
  
  # create columns for PODs
  info$POD_uni_centile <- NA
  info$POD_uni_mode <- NA
  info$POD_uni_lcrd <- NA
  info$POD_uni_gs <- NA
  
  info$POD_multi_gmd <- NA
  info$POD_multi_gmd_l <- NA
  info$POD_multi_gmd_u <- NA
  
  info$POD_multi_lmd <- NA
  info$POD_multi_lmd_l <- NA
  info$POD_multi_lmd_u <- NA
  
  info$POD_multi_pls <- NA
  info$POD_multi_pls_l <- NA
  info$POD_multi_pls_u <- NA
  
  for(k in c(1:length(analysisIDs))) {
    print(analysisIDs[k])
    
    bmds.temp <- bmds[bmds$analysisID == analysisIDs[k], ]
    if(grepl("httr", analysisIDs[k])) {
      platform <- "httr"
    } else {
      platform <- "htpp"
    }
  
    ## univariate PODs
    
    # compute 10th centile
    info$POD_uni_centile[k] <- centilePOD(bmds.temp$bmd, 0.1)
    
    # compute mode
    info$POD_uni_mode[k] <- modePOD(bmds.temp$bmd)
    
    # compute lcrd
    info$POD_uni_lcrd[k] <- lcrdPOD(bmds.temp$bmd)
    
    # compute pathway (go bp)
    if(platform == "httr") {
      gs.bmds <- bmds.temp$bmd
      names(gs.bmds) <- bmds.temp$gene.id
      info$POD_uni_gs[k] <- gsPOD(gs.list, universe, gs.bmds)
    } else {
      info$POD_uni_gs[k] <- NA
    }
    
    ## multivariate PODs
    
    # read in & process data matrix
    mat <- read_parquet(mat.paths[grep(analysisIDs[k], mat.paths)]) %>% as.data.frame()
    dose <- mat[1,-1] %>% unlist()
    rownames(mat) <- mat$probe_id
    mat <- mat[-1, -1] %>% t()
    
    # compute global mahalanobis (cov from all samples)
    if(platform == "httr"){
      gmd.rot <- httr.gmd$RotationMatrix
      gmd.cov <- httr.gmd$invCov
    } else {
      gmd.rot <- htpp.gmd$RotationMatrix
      gmd.cov <- htpp.gmd$invCov
    }
    
    gmd.rot <- gmd.rot[rownames(gmd.rot) %in% colnames(mat),]
    gmd.mat <- mat[,colnames(mat) %in% rownames(gmd.rot)]
    gmd.mat <- gmd.mat[,match(rownames(gmd.rot), colnames(gmd.mat))]
    
    gmd <- computeMahalanobisDistance(gmd.mat, gmd.rot, gmd.cov, as.character(dose), "0")
    gmd <- matrix(gmd, nrow = 1)
    rownames(gmd) <- "gmd"
    gmd.pod <- scoresPOD(gmd, dose)
    if(gmd.pod$all.pass){
      info$POD_multi_gmd[k] <- gmd.pod$bmd
      info$POD_multi_gmd_l[k] <- gmd.pod$bmdl
      info$POD_multi_gmd_u[k] <- gmd.pod$bmdu
    } else {
      info$POD_multi_gmd[k] <- NA
      info$POD_multi_gmd_l[k] <- NA
      info$POD_multi_gmd_u[k] <- NA
    }
    
    # compute local mahalanobis (cov from dataset only)
    lmd.cov <- prepMahalanobisDistances(mat, 0.95, as.character(dose))
    lmd <- computeMahalanobisDistance(mat, lmd.cov$RotationMatrix, lmd.cov$invCov, as.character(dose), "0")
    lmd <- matrix(lmd, nrow = 1)
    rownames(lmd) <- "lmd"
    lmd.pod <- scoresPOD(lmd, dose)
    if(lmd.pod$all.pass){
      info$POD_multi_lmd[k] <- lmd.pod$bmd
      info$POD_multi_lmd_l[k] <- lmd.pod$bmdl
      info$POD_multi_lmd_u[k] <- lmd.pod$bmdu
    } else {
      info$POD_multi_lmd[k] <- NA
      info$POD_multi_lmd_l[k] <- NA
      info$POD_multi_lmd_u[k] <- NA
    }
    
    # compute pls
    pls.comps <- computePLS(mat, as.character(dose), ncomp = 10)
    pls.bmds <- scoresPOD(pls.comps$scores, dose)
    pls.bmds <- pls.bmds[pls.bmds$all.pass, ]
    if(dim(pls.bmds)[1] > 0){
      pls.bmds$comp_num <- gsub("Comp ", "", pls.bmds$gene.id) %>% as.numeric()
      pls.bmds <- pls.bmds[order(pls.bmds$comp_num), ]
      
      info$POD_multi_pls[k] <- pls.bmds$bmd[1]
      info$POD_multi_pls_l[k] <- pls.bmds$bmdl[1]
      info$POD_multi_pls_u[k] <- pls.bmds$bmdu[1]
      
      # save loadings
      loadings <- pls.comps$loadings[,colnames(pls.comps$loadings) %in% pls.bmds$gene.id]
      if(platform == "httr"){
        httr.pls.loadings[[analysisIDs[k]]] <- loadings
      } else {
        htpp.pls.loadings[[analysisIDs[k]]] <- loadings
      }
      
    } else {
      info$POD_multi_pls[k] <- NA
      info$POD_multi_pls_l[k] <- NA
      info$POD_multi_pls_u[k] <- NA
    }
    
  }
  
  tpod_res <- rbind(tpod_res, info)
  
}


#### Compute differences in tPOD distributions across conditions ####

# Metrics: 
# KS p-value for tPOD dists between conditions, for each of the 11 chemicals (paired?)
# Paired difference of tPOD dist mean (1 test for all 11 chemicals) between conditions
# Paired difference of tPOD dist variance (1 test for all 11 chemicals) between conditions

# Across tPOD types for the same dataset

# Across datasets for the same tpod type

