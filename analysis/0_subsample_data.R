# Sub-sample matrices for all downstream analysis
# Jessica Ewald
# January 18, 2024


# The purpose of this script is to process and sub-sample both the transcriptomics 
# (httr) and Cell Painting (htpp) matrices to create many versions, each with fewer 
# and randomly sampled replicates. 


#### Set up ####

library(dplyr)
library(edgeR)
library(arrow)

# relative path to data submodule
data.path <- "../2024_01_18_Dose_response_methods-data/data/"

# read in files & convert to proper type
httr_key <- read_parquet(paste0(data.path, "1_orig_data/httr_key.parquet")) %>% as.data.frame()

httr_counts <- read_parquet(paste0(data.path, "1_orig_data/httr_counts.parquet")) %>% as.data.frame()
rownames(httr_counts) <- httr_counts$probe_id
httr_counts <- httr_counts[,-1]
httr_counts <- as.matrix(httr_counts)

htpp_chem <- read_parquet(paste0(data.path, "1_orig_data/htpp_chem.parquet")) %>% as.data.frame()

htpp_feature <- read_parquet(paste0(data.path, "1_orig_data/htpp_feature.parquet")) %>% as.data.frame()

htpp_global_mah <- read_parquet(paste0(data.path, "1_orig_data/htpp_global_mah.parquet")) %>% as.data.frame()

htpp_well_norm <- read_parquet(paste0(data.path, "1_orig_data/htpp_well_norm.parquet")) %>% as.data.frame()


#### Process transcriptomics data ####

# first, retain only samples with QA/QC = OK
httr_counts <- httr_counts[,colnames(httr_counts) %in% httr_key$sample_id[httr_key$qc_flag == "OK"]]

# next, separate matrix into a separate matrix for each chemical (filter, convert to logCPM)
chems <- unique(httr_key$chem_name[!is.na(httr_key$chem_name)])
httr_lcpm <- list()

for(i in c(1:length(chems))){
  # grab ctrls + samples from specific chemical
  temp.samps <- httr_key$sample_id[httr_key$chem_name == chems[i] | is.na(httr_key$chem_name)]
  temp.count <- httr_counts[,colnames(httr_counts) %in% temp.samps]
  
  # filter based on abundance
  keep.feat <- apply(temp.count, 1, mean) > 5
  temp.count <- temp.count[keep.feat, ]
  
  # convert to log-cpm
  nf <- calcNormFactors(temp.count, method="RLE")
  temp.lcpm <- cpm(temp.count, lib.size=colSums(temp.count)*nf, log = TRUE, prior.count = 1)
  
  # re-order by dose and save dose
  temp.meta <- httr_key[httr_key$sample_id %in% colnames(temp.lcpm),]
  temp.meta <- temp.meta[order(temp.meta$conc), ]
  temp.conc <- temp.meta$conc
  names(temp.conc) <- temp.meta$sample_id
  temp.lcpm <- temp.lcpm[,match(temp.meta$sample_id, colnames(temp.lcpm))]
  identical(colnames(temp.lcpm), temp.meta$sample_id) # double check in order
  
  httr_lcpm[[i]] <- list(lcpm = temp.lcpm, conc = temp.conc)
}
names(httr_lcpm) <- chems


#### Sub-sample transcriptomics data ####

# next, we need to sub-sample to create multiple matrices per chemical and per replicate number
# sample with replacement when less than 3 per group after removing QC fails, keep all ctrls

for(i in c(1:length(chems))){
  chem <- chems[i]
  conc <- httr_lcpm[[chem]]$conc
  lcpm <- httr_lcpm[[chem]]$lcpm
  
  ctrl <- lcpm[,colnames(lcpm) %in% names(conc)[conc == 0]]
  concs <- conc[conc > 0] # remove ctrls
  conc.vals <- unique(concs)
  
  # create for 2 reps
  for(k in c(1:30)){
    sel.samps <- lapply(conc.vals, function(x){
      sample(names(conc)[conc == x], size = 2, replace = FALSE)
    }) %>% unlist()
    
    temp.lcpm <- cbind(ctrl, lcpm[,colnames(lcpm) %in% sel.samps])
    temp.conc <- conc[names(conc) %in% colnames(temp.lcpm)]
    
    # write out as parquet
    temp.lcpm <- rbind(temp.conc, temp.lcpm)
    rownames(temp.lcpm)[1] <- "Dose"
    temp.lcpm <- cbind(rownames(temp.lcpm), temp.lcpm)
    colnames(temp.lcpm)[1] <- "probe_id"
    write_parquet(temp.lcpm, 
                  paste0(data.path, "2_sampled_data/httr/2_reps/", chem, "_2reps_", k, ".parquet"))
    
    # write out in BMDExpress format for WTT
    write.table(temp.lcpm, paste0(data.path, "2b_BMDExpress_input/", chem, "_2reps_", k, ".txt"), sep="\t", col.names = NA)
  }
  
  # create for 3 reps
  for(k in c(1:30)){
    sel.samps <- lapply(conc.vals, function(x){
      if(length(names(conc)[conc == x]) > 2){
        sample(names(conc)[conc == x], size = 3, replace = FALSE)
      } else {
        sample(names(conc)[conc == x], size = 3, replace = TRUE)
      }
    }) %>% unlist()
    
    temp.lcpm <- cbind(ctrl, lcpm[,colnames(lcpm) %in% sel.samps])
    temp.conc <- conc[names(conc) %in% colnames(temp.lcpm)]
    
    # write out as parquet
    temp.lcpm <- rbind(temp.conc, temp.lcpm)
    rownames(temp.lcpm)[1] <- "Dose"
    temp.lcpm <- cbind(rownames(temp.lcpm), temp.lcpm)
    colnames(temp.lcpm)[1] <- "probe_id"
    write_parquet(temp.lcpm, 
                  paste0(data.path, "2_sampled_data/httr/3_reps/", chem, "_3reps_", k, ".parquet"))
    
    # write out in BMDExpress format for WTT
    write.table(temp.lcpm, paste0(data.path, "2b_BMDExpress_input/", chem, "_3reps_", k, ".txt"), sep="\t", col.names = NA)
  }
}


#### Process Cell Painting data ####
# next, separate into a separate matrix for each chemical
htpp_well_norm$conc[is.na(htpp_well_norm$conc)] <- 0
htpp_mats <- list()

for(i in c(1:length(chems))){
  
  # grab ctrls + samples from specific chemical
  temp.samps <- htpp_well_norm$sample_id[htpp_well_norm$chem_id == chems[i] | htpp_well_norm$chem_id == "Dimethyl Sulfoxide"]
  temp.mat <- htpp_well_norm[htpp_well_norm$sample_id %in% temp.samps, ]
  rownames(temp.mat) <- temp.mat$sample_id
  temp.conc <- temp.mat$conc
  temp.mat <- temp.mat[,-c(1:11)]
  temp.mat <- t(temp.mat) %>% as.data.frame()
  
  # re-order by dose and save dose
  names(temp.conc) <- colnames(temp.mat)
  temp.conc <- sort(temp.conc)
  temp.mat <- temp.mat[,match(names(temp.conc), colnames(temp.mat))]
  identical(colnames(temp.mat), names(temp.conc)) # double check in order
  
  htpp_mats[[i]] <- list(dat = temp.mat, conc = temp.conc)
}
names(htpp_mats) <- chems


#### Sub-sample Cell Painting data ####

for(i in c(1:length(chems))){
  chem <- chems[i]
  conc <- htpp_mats[[chem]]$conc
  dat <- htpp_mats[[chem]]$dat
  
  ctrl <- dat[,colnames(dat) %in% names(conc)[conc == 0]]
  concs <- conc[conc > 0] # remove ctrls
  conc.vals <- unique(concs)
  
  # create for 2 reps
  for(k in c(1:30)){
    sel.samps <- lapply(conc.vals, function(x){
      sample(names(conc)[conc == x], size = 2, replace = FALSE)
    }) %>% unlist()
    
    temp.dat <- cbind(ctrl, dat[,colnames(dat) %in% sel.samps])
    temp.conc <- conc[names(conc) %in% colnames(temp.dat)]
    
    # write out as parquet
    temp.dat <- rbind(temp.conc, temp.dat)
    rownames(temp.dat)[1] <- "Dose"
    temp.dat <- cbind(rownames(temp.dat), temp.dat)
    colnames(temp.dat)[1] <- "probe_id"
    write_parquet(temp.dat, 
                  paste0(data.path, "2_sampled_data/htpp/2_reps/", chem, "_2reps_", k, "_htpp.parquet"))
    
  }
  
  # create for 3 reps
  for(k in c(1:30)){
    sel.samps <- lapply(conc.vals, function(x){
      if(length(names(conc)[conc == x]) > 2){
        sample(names(conc)[conc == x], size = 3, replace = FALSE)
      } else {
        sample(names(conc)[conc == x], size = 3, replace = TRUE)
      }
    }) %>% unlist()
    
    temp.dat <- cbind(ctrl, dat[,colnames(dat) %in% sel.samps])
    temp.conc <- conc[names(conc) %in% colnames(temp.dat)]
    
    # write out as parquet
    temp.dat <- rbind(temp.conc, temp.dat)
    rownames(temp.dat)[1] <- "Dose"
    temp.dat <- cbind(rownames(temp.dat), temp.dat)
    colnames(temp.dat)[1] <- "probe_id"
    write_parquet(temp.dat, 
                  paste0(data.path, "2_sampled_data/htpp/3_reps/", chem, "_3reps_", k, "_htpp.parquet"))
  }
}

