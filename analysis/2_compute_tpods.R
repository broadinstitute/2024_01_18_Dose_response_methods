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
               paste0(data.path, "4_bmd_results/htpp_bmd_wtt.parquet"),
               paste0(data.path, "4_bmd_results/httr_fc/httr_bmd_fc1.parquet"),
               paste0(data.path, "4_bmd_results/httr_fc/httr_bmd_fc2.parquet"),
               paste0(data.path, "4_bmd_results/httr_fc/httr_bmd_fc3.parquet"),
               paste0(data.path, "4_bmd_results/httr_fc/httr_bmd_fc4.parquet"),
               paste0(data.path, "4_bmd_results/httr_fc/httr_bmd_fc5.parquet"))

names(bmd.paths) <- c("httr_aov", "httr_wtt", "httr_s1500", "httr_nomic", "htpp_all", "htpp_aov", "htpp_wtt", 
                      "httr_fc1", "httr_fc2", "httr_fc3", "httr_fc4", "httr_fc5")
feature.set <- c("ANOVA", "WTT", "S1500", "Nomic", "All", "ANOVA", "WTT", "FC1", "FC2", "FC3", "FC4", "FC5")

# get gene set info
gs.list <- readRDS("../2024_01_18_Dose_response_methods-data/data/gs_libraries/go_bp.rds")[["sets"]]
probe.map <- read.csv("../2024_01_18_Dose_response_methods-data/data/probemaps/100739_homo_sapiens_wt_probemap.csv")
universe.httr <- probe.map$PROBE_NAME %>% unique()

s1500.map <- read.csv("../2024_01_18_Dose_response_methods-data/data/probemaps/100766_homo_sapiens_s1500_probemap.csv")
universe.s1500 <- s1500.map$PROBE_NAME

nomic.map <- read.csv("../2024_01_18_Dose_response_methods-data/data/probemaps/Nomic_probes_S1500_jan2024.csv")
nomic.map <- nomic.map[nomic.map$Status != "On roadmap", ]
nomic.map <- merge(nomic.map, probe.map, by = "ENTREZ_ID", all = FALSE)
universe.nomic <- unique(nomic.map$PROBE_NAME)

universe.htpp <- read_parquet("../2024_01_18_Dose_response_methods-data/data/2_sampled_data/htpp/3_reps/Actinomycin D_3reps_1_htpp.parquet")$probe_id[-1]

#### Compute univariate tPODs ####

tpod.res <- data.frame()
for(i in c(5:length(bmd.paths))){
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
  info$POD_uni_rs <- NA
  
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
    mode <- modePOD(log10(bmds.temp$bmd))
    if(!is.na(mode)){
      info$POD_uni_mode[k] <- 10^mode
    } else {
      info$POD_uni_mode[k] <- NA
    }
    
    # compute lcrd
    info$POD_uni_lcrd[k] <- lcrdPOD(bmds.temp$bmd)
    
    # compute gene set and random set POD
    gs.bmds <- bmds.temp$bmd
    names(gs.bmds) <- bmds.temp$gene.id
    
    if(feature.set[i] == "S1500"){
      info$POD_uni_rs[k] <- rsPOD(universe.s1500, gs.bmds, 100)[1]
      info$POD_uni_gs[k] <- gsPOD(gs.list, universe.s1500, gs.bmds, 3, 0.05)
    } else if (feature.set[i] == "Nomic"){
      info$POD_uni_rs[k] <- rsPOD(universe.nomic, gs.bmds, 100)[1]
      info$POD_uni_gs[k] <- gsPOD(gs.list, universe.nomic, gs.bmds, 3, 0.05)
    } else if(platform == "httr") {
      info$POD_uni_rs[k] <- rsPOD(universe.httr, gs.bmds, 100)[1]
      info$POD_uni_gs[k] <- gsPOD(gs.list, universe.httr, gs.bmds, 3, 0.05)
    } else {
      info$POD_uni_rs[k] <- rsPOD(universe.htpp, gs.bmds, 100)[1]
      info$POD_uni_gs[k] <- NA
    }
  }
  
  tpod.res <- rbind(tpod.res, info)
  
}
write_parquet(tpod.res, paste0(data.path, "5_tpod_results/univar_tpods.parquet"))


# HTTR ANOVA, HTTR WTT, HTPP All, HTPP ANOVA, HTPP WTT all fine (except for NAs - probably just a few NAs)
# Don't need Nomic
# Re-do S1500
tpod <- read_parquet("/Users/jessicaewald/NetbeansProjects/2024_01_18_Dose_response_methods/2024_01_18_Dose_response_methods-data/data/5_tpod_results/univar_tpods.parquet") %>% as.data.frame()

### wipe ones that didn't work
tpod <- tpod[tpod$Features != "Nomic", ]
tpod$POD_uni_gs[tpod$Features == "S1500"] <- NA
tpod$POD_uni_rs[tpod$Features == "S1500"] <- NA

## re-do - collect bmd.res into single df
bmd.res <- read_parquet(bmd.paths["httr_aov"]) %>% as.data.frame()
bmd.res$Features <- "ANOVA"

bmd2 <- read_parquet(bmd.paths["httr_wtt"]) %>% as.data.frame()
bmd2$Features <- "WTT"
bmd.res <- rbind(bmd.res, bmd2)

bmd2 <- read_parquet(bmd.paths["httr_s1500"]) %>% as.data.frame()
bmd2$Features <- "S1500"
bmd.res <- rbind(bmd.res, bmd2)

bmd2 <- read_parquet(bmd.paths["htpp_all"]) %>% as.data.frame()
bmd2$Features <- "All"
bmd.res <- rbind(bmd.res, bmd2)

bmd2 <- read_parquet(bmd.paths["htpp_aov"]) %>% as.data.frame()
bmd2$Features <- "ANOVA"
bmd.res <- rbind(bmd.res, bmd2)

bmd2 <- read_parquet(bmd.paths["htpp_wtt"]) %>% as.data.frame()
bmd2$Features <- "WTT"
bmd.res <- rbind(bmd.res, bmd2)

# go through tpod and fill in NAs for gs and rs PODs
for(i in c(1:dim(tpod)[1])){
  
  if(tpod$Features[i] == "S1500"){
    print(i)
    bmds.temp <- bmd.res[bmd.res$analysisID == tpod$AnalysisID[i] & bmd.res$Features == "S1500", ]
    gs.bmds <- bmds.temp$bmd
    names(gs.bmds) <- bmds.temp$gene.id
    
    tpod$POD_uni_rs[i] <- rsPOD(universe.s1500, gs.bmds, 100)[1]
    tpod$POD_uni_gs[i] <- gsPOD(gs.list, universe.s1500, gs.bmds, 3, 0.05)
  } else if(is.na(tpod$POD_uni_rs[i])){
    print(i)
    bmds.temp <- bmd.res[bmd.res$analysisID == tpod$AnalysisID[i] & bmd.res$Features == tpod$Features[i], ]
    gs.bmds <- bmds.temp$bmd
    names(gs.bmds) <- bmds.temp$gene.id
    
    if(tpod$Platform[i] == "httr") {
      tpod$POD_uni_rs[i] <- rsPOD(universe.httr, gs.bmds, 100)[1]
    } else {
      tpod$POD_uni_rs[i] <- rsPOD(universe.htpp, gs.bmds, 100)[1]
    }
  }
}
write_parquet(tpod, paste0(data.path, "5_tpod_results/univar_tpods.parquet"))

#### Compute multivariate PODs ####

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

# file paths
mat.paths <- c(list.files(paste0(data.path, "2_sampled_data/httr/2_reps"), full.names = T),
               list.files(paste0(data.path, "2_sampled_data/httr/3_reps"), full.names = T),
               list.files(paste0(data.path, "2_sampled_data/htpp/2_reps"), full.names = T),
               list.files(paste0(data.path, "2_sampled_data/htpp/3_reps"), full.names = T))
mat.names <- c(list.files(paste0(data.path, "2_sampled_data/httr/2_reps"), full.names = F),
               list.files(paste0(data.path, "2_sampled_data/httr/3_reps"), full.names = F),
               list.files(paste0(data.path, "2_sampled_data/htpp/2_reps"), full.names = F),
               list.files(paste0(data.path, "2_sampled_data/htpp/3_reps"), full.names = F)) %>% gsub("\\.parquet", "", .)

tpods.multi <- data.frame(AnalysisID = mat.names,
                          POD_multi_gmd = NA,
                          POD_multi_gmd_l = NA,
                          POD_multi_gmd_u = NA,
                          POD_multi_lmd = NA,
                          POD_multi_lmd_l = NA,
                          POD_multi_lmd_u = NA,
                          POD_multi_pls = NA,
                          POD_multi_pls_l = NA,
                          POD_multi_pls_u = NA)
pls.loadings <- list()

for(i in c(1:length(mat.paths))){
  ## multivariate PODs
  print(mat.names[i])
  
  # set platform
  if(grepl("httr", mat.names[i])) {
    platform <- "httr"
  } else {
    platform <- "htpp"
  }
  
  # read in & process data matrix
  mat <- read_parquet(mat.paths[i]) %>% as.data.frame()
  dose <- mat[1,-1] %>% unlist()
  rownames(mat) <- mat$probe_id
  mat <- mat[-1, -1] %>% t()
  
  # compute global mahalanobis (cov from all samples)
  if(platform == "httr"){
    gmd.rot <- httr.gmd$RotationMatrix
    gmd.cov <- httr.gmd$invCov
  } else {
    mat <- mat[,!(colnames(mat) %in% c("f_1073", "f_1074"))]
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
    tpods.multi$POD_multi_gmd[i] <- gmd.pod$bmd
    tpods.multi$POD_multi_gmd_l[i] <- gmd.pod$bmdl
    tpods.multi$POD_multi_gmd_u[i] <- gmd.pod$bmdu
  } else {
    tpods.multi$POD_multi_gmd[i] <- NA
    tpods.multi$POD_multi_gmd_l[i] <- NA
    tpods.multi$POD_multi_gmd_u[i] <- NA
  }

  # compute local mahalanobis (cov from dataset only)
  lmd.cov <- prepMahalanobisDistances(mat, 0.95, as.character(dose))
  lmd <- computeMahalanobisDistance(mat, lmd.cov$RotationMatrix, lmd.cov$invCov, as.character(dose), "0")
  lmd <- matrix(lmd, nrow = 1)
  rownames(lmd) <- "lmd"
  lmd.pod <- scoresPOD(lmd, dose)
  if(lmd.pod$all.pass){
    tpods.multi$POD_multi_lmd[i] <- lmd.pod$bmd
    tpods.multi$POD_multi_lmd_l[i] <- lmd.pod$bmdl
    tpods.multi$POD_multi_lmd_u[i] <- lmd.pod$bmdu
  } else {
    tpods.multi$POD_multi_lmd[i] <- NA
    tpods.multi$POD_multi_lmd_l[i] <- NA
    tpods.multi$POD_multi_lmd_u[i] <- NA
  }
  
  # compute pls
  pls.comps <- computePLS(mat, as.character(dose), ncomp = 10)
  pls.bmds <- scoresPOD(pls.comps$scores, dose)
  pls.bmds <- pls.bmds[pls.bmds$all.pass, ]
  if(dim(pls.bmds)[1] > 0){
    pls.bmds$comp_num <- gsub("Comp ", "", pls.bmds$gene.id) %>% as.numeric()
    pls.bmds <- pls.bmds[order(pls.bmds$comp_num), ]
    
    tpods.multi$POD_multi_pls[i] <- pls.bmds$bmd[1]
    tpods.multi$POD_multi_pls_l[i] <- pls.bmds$bmdl[1]
    tpods.multi$POD_multi_pls_u[i] <- pls.bmds$bmdu[1]
    
    # save loadings
    loadings <- pls.comps$loadings[,colnames(pls.comps$loadings) %in% pls.bmds$gene.id]
    pls.loadings[[mat.names[i]]] <- loadings
    
  } else {
    tpods.multi$POD_multi_pls[i] <- NA
    tpods.multi$POD_multi_pls_l[i] <- NA
    tpods.multi$POD_multi_pls_u[i] <- NA
  }
}

#write_parquet(tpods.multi, paste0(data.path, "5_tpod_results/multivar_tpods.parquet"))
saveRDS(pls.loadings, paste0(data.path, "pls_loadings/pls_loadings.rds"))

#### Compute multivariate PODs for reduced transcriptomes ####
s1500 <- read.csv("/Users/jessicaewald/NetbeansProjects/2024_01_18_Dose_response_methods/2024_01_18_Dose_response_methods-data/data/probemaps/100766_homo_sapiens_s1500_probemap.csv")
nomic <- read.csv("/Users/jessicaewald/NetbeansProjects/2024_01_18_Dose_response_methods/2024_01_18_Dose_response_methods-data/data/probemaps/Nomic_probes_S1500_jan2024.csv")

nomic <- merge(nomic, s1500[,c("ENTREZ_ID", "PROBE_NAME")], by = "ENTREZ_ID")

# get global MD pre-reqs
# full httr matrix
httr <- read_parquet("../2024_01_18_Dose_response_methods-data/data/1_orig_data/httr_counts.parquet") %>% as.data.frame()
httr.meta <- read_parquet("../2024_01_18_Dose_response_methods-data/data/1_orig_data/httr_key.parquet") %>% as.data.frame()
rownames(httr) <- httr$probe_id
httr <- httr[,-1]
httr <- as.matrix(httr)
httr <- httr[,colnames(httr) %in% httr.meta$sample_id[httr.meta$qc_flag == "OK"]]
httr.meta <- httr.meta[httr.meta$qc_flag == "OK", ]

# process s1500
s1500.httr <- httr[rownames(httr) %in% s1500$PROBE_NAME, ]
keep.feat <- apply(s1500.httr, 1, mean) > 5
s1500.httr <- s1500.httr[keep.feat, ]
nf <- calcNormFactors(s1500.httr, method="RLE")
s1500.httr <- cpm(s1500.httr, lib.size=colSums(s1500.httr)*nf, log = TRUE, prior.count = 1)
s1500.httr <- t(s1500.httr)

# process nomic
nomic.httr <- httr[rownames(httr) %in% nomic$PROBE_NAME, ]
keep.feat <- apply(nomic.httr, 1, mean) > 5
nomic.httr <- nomic.httr[keep.feat, ]
nf <- calcNormFactors(nomic.httr, method="RLE")
nomic.httr <- cpm(nomic.httr, lib.size=colSums(nomic.httr)*nf, log = TRUE, prior.count = 1)
nomic.httr <- t(nomic.httr)

# compute MD stats
s1500.gmd <- prepMahalanobisDistances(s1500.httr, 0.95, paste0(httr.meta$chem_name, "_", httr.meta$dose_level))
nomic.gmd <- prepMahalanobisDistances(nomic.httr, 0.95, paste0(httr.meta$chem_name, "_", httr.meta$dose_level))


# file paths
mat.paths <- c(list.files(paste0(data.path, "2_sampled_data/httr/2_reps"), full.names = T),
               list.files(paste0(data.path, "2_sampled_data/httr/3_reps"), full.names = T))
mat.names <- c(list.files(paste0(data.path, "2_sampled_data/httr/2_reps"), full.names = F),
               list.files(paste0(data.path, "2_sampled_data/httr/3_reps"), full.names = F)) %>% gsub("\\.parquet", "", .)

tpods.multi.s1500 <- data.frame(AnalysisID = mat.names,
                                POD_multi_gmd = NA,
                                POD_multi_gmd_l = NA,
                                POD_multi_gmd_u = NA,
                                POD_multi_lmd = NA,
                                POD_multi_lmd_l = NA,
                                POD_multi_lmd_u = NA,
                                POD_multi_pls = NA,
                                POD_multi_pls_l = NA,
                                POD_multi_pls_u = NA)

tpods.multi.nomic <- data.frame(AnalysisID = mat.names,
                                POD_multi_gmd = NA,
                                POD_multi_gmd_l = NA,
                                POD_multi_gmd_u = NA,
                                POD_multi_lmd = NA,
                                POD_multi_lmd_l = NA,
                                POD_multi_lmd_u = NA,
                                POD_multi_pls = NA,
                                POD_multi_pls_l = NA,
                                POD_multi_pls_u = NA)

pls.loadings.s1500 <- list()
pls.loadings.nomic <- list()

##### S1500
for(i in c(1:length(mat.paths))){
  
  ## multivariate PODs
  print(mat.names[i])
  
  # read in & process data matrix
  mat <- read_parquet(mat.paths[i]) %>% as.data.frame()
  dose <- mat[1,-1] %>% unlist()
  rownames(mat) <- mat$probe_id
  mat <- mat[-1, -1] %>% t()
  mat <- mat[,colnames(mat) %in% s1500$PROBE_NAME]
  
  # compute global mahalanobis (cov from all samples)
  gmd.rot <- s1500.gmd$RotationMatrix
  gmd.cov <- s1500.gmd$invCov
  
  gmd.rot <- gmd.rot[rownames(gmd.rot) %in% colnames(mat),]
  gmd.mat <- mat[,colnames(mat) %in% rownames(gmd.rot)]
  gmd.mat <- gmd.mat[,match(rownames(gmd.rot), colnames(gmd.mat))]
  
  gmd <- computeMahalanobisDistance(gmd.mat, gmd.rot, gmd.cov, as.character(dose), "0")
  gmd <- matrix(gmd, nrow = 1)
  rownames(gmd) <- "gmd"
  gmd.pod <- scoresPOD(gmd, dose)
  if(gmd.pod$all.pass){
    tpods.multi.s1500$POD_multi_gmd[i] <- gmd.pod$bmd
    tpods.multi.s1500$POD_multi_gmd_l[i] <- gmd.pod$bmdl
    tpods.multi.s1500$POD_multi_gmd_u[i] <- gmd.pod$bmdu
  } else {
    tpods.multi.s1500$POD_multi_gmd[i] <- NA
    tpods.multi.s1500$POD_multi_gmd_l[i] <- NA
    tpods.multi.s1500$POD_multi_gmd_u[i] <- NA
  }
  
  # compute local mahalanobis (cov from dataset only)
  lmd.cov <- prepMahalanobisDistances(mat, 0.95, as.character(dose))
  lmd <- computeMahalanobisDistance(mat, lmd.cov$RotationMatrix, lmd.cov$invCov, as.character(dose), "0")
  lmd <- matrix(lmd, nrow = 1)
  rownames(lmd) <- "lmd"
  lmd.pod <- scoresPOD(lmd, dose)
  if(lmd.pod$all.pass){
    tpods.multi.s1500$POD_multi_lmd[i] <- lmd.pod$bmd
    tpods.multi.s1500$POD_multi_lmd_l[i] <- lmd.pod$bmdl
    tpods.multi.s1500$POD_multi_lmd_u[i] <- lmd.pod$bmdu
  } else {
    tpods.multi.s1500$POD_multi_lmd[i] <- NA
    tpods.multi.s1500$POD_multi_lmd_l[i] <- NA
    tpods.multi.s1500$POD_multi_lmd_u[i] <- NA
  }
  
  # compute pls
  pls.comps <- computePLS(mat, as.character(dose), ncomp = 10)
  pls.bmds <- scoresPOD(pls.comps$scores, dose)
  pls.bmds <- pls.bmds[pls.bmds$all.pass, ]
  if(dim(pls.bmds)[1] > 0){
    pls.bmds$comp_num <- gsub("Comp ", "", pls.bmds$gene.id) %>% as.numeric()
    pls.bmds <- pls.bmds[order(pls.bmds$comp_num), ]
    
    tpods.multi.s1500$POD_multi_pls[i] <- pls.bmds$bmd[1]
    tpods.multi.s1500$POD_multi_pls_l[i] <- pls.bmds$bmdl[1]
    tpods.multi.s1500$POD_multi_pls_u[i] <- pls.bmds$bmdu[1]
    
    # save loadings
    loadings <- pls.comps$loadings[,colnames(pls.comps$loadings) %in% pls.bmds$gene.id]
    pls.loadings.s1500[[mat.names[i]]] <- loadings
    
  } else {
    tpods.multi.s1500$POD_multi_pls[i] <- NA
    tpods.multi.s1500$POD_multi_pls_l[i] <- NA
    tpods.multi.s1500$POD_multi_pls_u[i] <- NA
  }
}
write_parquet(tpods.multi.s1500, paste0(data.path, "5_tpod_results/multivar_tpods_s1500.parquet"))
saveRDS(pls.loadings.s1500, paste0(data.path, "pls_loadings/pls_loadings_s1500.rds"))

##### NOMIC
for(i in c(1:length(mat.paths))){
  
  ## multivariate PODs
  print(mat.names[i])
  
  # read in & process data matrix
  mat <- read_parquet(mat.paths[i]) %>% as.data.frame()
  dose <- mat[1,-1] %>% unlist()
  rownames(mat) <- mat$probe_id
  mat <- mat[-1, -1] %>% t()
  mat <- mat[,colnames(mat) %in% nomic$PROBE_NAME]
  
  # compute global mahalanobis (cov from all samples)
  gmd.rot <- nomic.gmd$RotationMatrix
  gmd.cov <- nomic.gmd$invCov
  
  gmd.rot <- gmd.rot[rownames(gmd.rot) %in% colnames(mat),]
  gmd.mat <- mat[,colnames(mat) %in% rownames(gmd.rot)]
  gmd.mat <- gmd.mat[,match(rownames(gmd.rot), colnames(gmd.mat))]
  
  gmd <- computeMahalanobisDistance(gmd.mat, gmd.rot, gmd.cov, as.character(dose), "0")
  gmd <- matrix(gmd, nrow = 1)
  rownames(gmd) <- "gmd"
  gmd.pod <- scoresPOD(gmd, dose)
  if(gmd.pod$all.pass){
    tpods.multi.nomic$POD_multi_gmd[i] <- gmd.pod$bmd
    tpods.multi.nomic$POD_multi_gmd_l[i] <- gmd.pod$bmdl
    tpods.multi.nomic$POD_multi_gmd_u[i] <- gmd.pod$bmdu
  } else {
    tpods.multi.nomic$POD_multi_gmd[i] <- NA
    tpods.multi.nomic$POD_multi_gmd_l[i] <- NA
    tpods.multi.nomic$POD_multi_gmd_u[i] <- NA
  }
  
  # compute local mahalanobis (cov from dataset only)
  lmd.cov <- prepMahalanobisDistances(mat, 0.95, as.character(dose))
  lmd <- computeMahalanobisDistance(mat, lmd.cov$RotationMatrix, lmd.cov$invCov, as.character(dose), "0")
  lmd <- matrix(lmd, nrow = 1)
  rownames(lmd) <- "lmd"
  lmd.pod <- scoresPOD(lmd, dose)
  if(lmd.pod$all.pass){
    tpods.multi.nomic$POD_multi_lmd[i] <- lmd.pod$bmd
    tpods.multi.nomic$POD_multi_lmd_l[i] <- lmd.pod$bmdl
    tpods.multi.nomic$POD_multi_lmd_u[i] <- lmd.pod$bmdu
  } else {
    tpods.multi.nomic$POD_multi_lmd[i] <- NA
    tpods.multi.nomic$POD_multi_lmd_l[i] <- NA
    tpods.multi.nomic$POD_multi_lmd_u[i] <- NA
  }
  
  # compute pls
  pls.comps <- computePLS(mat, as.character(dose), ncomp = 10)
  pls.bmds <- scoresPOD(pls.comps$scores, dose)
  pls.bmds <- pls.bmds[pls.bmds$all.pass, ]
  if(dim(pls.bmds)[1] > 0){
    pls.bmds$comp_num <- gsub("Comp ", "", pls.bmds$gene.id) %>% as.numeric()
    pls.bmds <- pls.bmds[order(pls.bmds$comp_num), ]
    
    tpods.multi.nomic$POD_multi_pls[i] <- pls.bmds$bmd[1]
    tpods.multi.nomic$POD_multi_pls_l[i] <- pls.bmds$bmdl[1]
    tpods.multi.nomic$POD_multi_pls_u[i] <- pls.bmds$bmdu[1]
    
    # save loadings
    loadings <- pls.comps$loadings[,colnames(pls.comps$loadings) %in% pls.bmds$gene.id]
    pls.loadings.nomic[[mat.names[i]]] <- loadings
    
  } else {
    tpods.multi.nomic$POD_multi_pls[i] <- NA
    tpods.multi.nomic$POD_multi_pls_l[i] <- NA
    tpods.multi.nomic$POD_multi_pls_u[i] <- NA
  }
}

write_parquet(tpods.multi.nomic, paste0(data.path, "5_tpod_results/multivar_tpods_nomic.parquet"))
saveRDS(pls.loadings.nomic, paste0(data.path, "pls_loadings/pls_loadings_nomic.rds"))


