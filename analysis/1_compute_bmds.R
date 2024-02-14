# Compute BMDs for all genes of interest
# January 19, 2024
# Jessica Ewald

# The purpose of this script is to compute all the feature BMDs that will be used 
# in the rest of the analysis. These include ANOVA, WTT, S1500, and Nomic probes 
# for the transcriptomics data, and all Cell Painting features. 

#### Setup ####

library(dplyr)
library(arrow)
source("../fastbmdR/fastbmdR_main.R")
source("../fastbmdR/fastbmdR_utils.R")

# relative path to data submodule
data.path <- "../2024_01_18_Dose_response_methods-data/data/"

# get lists of paths and analysis IDs 
httr.paths <- c(list.files(paste0(data.path, "2_sampled_data/httr/2_reps/"), full.names = TRUE),
                list.files(paste0(data.path, "2_sampled_data/httr/3_reps/"), full.names = TRUE))
httr.nms <- c(list.files(paste0(data.path, "2_sampled_data/httr/2_reps/"), full.names = FALSE),
              list.files(paste0(data.path, "2_sampled_data/httr/3_reps/"), full.names = FALSE))
httr.nms <- gsub(".parquet", "", httr.nms)

htpp.paths <- c(list.files(paste0(data.path, "2_sampled_data/htpp/2_reps/"), full.names = TRUE),
                list.files(paste0(data.path, "2_sampled_data/htpp/3_reps/"), full.names = TRUE))
htpp.nms <- c(list.files(paste0(data.path, "2_sampled_data/htpp/2_reps/"), full.names = FALSE),
              list.files(paste0(data.path, "2_sampled_data/htpp/3_reps/"), full.names = FALSE))
htpp.nms <- gsub(".parquet", "", htpp.nms)


#### Get feature filters #####

## Read in WTT results from BMDExpress
process_wtt <- function(wtt.path){
  wtt <- read.table(wtt.path, sep="\t", fill = TRUE, header = TRUE)
  colnames(wtt) <- c("analysis", "probe", "pval", "fdr", "max_fc")
  wtt$analysisID <- gsub("_williams.*", "", wtt$analysis)
  return(wtt)
}
httr.wtt.filter.2 <- process_wtt(paste0(data.path, "2c_WTT_results/processed/httr_2reps_wtt.txt"))
httr.wtt.filter.3 <- process_wtt(paste0(data.path, "2c_WTT_results/processed/httr_3reps_wtt.txt"))
httr.wtt.filter <- rbind(httr.wtt.filter.2, httr.wtt.filter.3)

htpp.wtt.filter.2 <- process_wtt(paste0(data.path, "2c_WTT_results/processed/htpp_2reps_wtt.txt"))
htpp.wtt.filter.3 <- process_wtt(paste0(data.path, "2c_WTT_results/processed/htpp_3reps_wtt.txt"))
htpp.wtt.filter <- rbind(htpp.wtt.filter.2, htpp.wtt.filter.3)

# write out results to parquet
write_parquet(httr.wtt.filter, paste0(data.path, "3_filtered_features/httr_wtt_filter.parquet"))
write_parquet(htpp.wtt.filter, paste0(data.path, "3_filtered_features/htpp_wtt_filter.parquet"))

## Create ANOVA filter lists for httr and htpp
httr.aov.filter <- data.frame()
for(i in c(1:length(httr.paths))){
  analysisID <- httr.nms[i]
  print(analysisID)
  
  # read in and process httr data
  lcpm <- read_parquet(httr.paths[i]) %>% as.data.frame()
  conc <- lcpm[1,-1] %>% unlist()
  rownames(lcpm) <- lcpm$probe_id
  lcpm <- lcpm[-1,-1]
  
  # perform prefiltering - simple filter (ANOVA)
  pvals <- apply(lcpm, 1, function(x){res <- anova(aov(x ~ conc))$`Pr(>F)`[1]})
  aov.df <- data.frame(probe = names(pvals), pval = unname(pvals), analysisID = analysisID)
  aov.df <- aov.df[aov.df$pval < 0.05, ]
  httr.aov.filter <- rbind(httr.aov.filter, aov.df)
}

htpp.aov.filter <- data.frame()
for(i in c(1:length(htpp.paths))){
  analysisID <- htpp.nms[i]
  print(analysisID)
  
  # read in and process httr data
  dat <- read_parquet(htpp.paths[i]) %>% as.data.frame()
  conc <- dat[1,-1] %>% unlist()
  rownames(dat) <- dat$probe_id
  dat <- dat[-1,-1]
  
  # perform prefiltering - simple filter (ANOVA)
  pvals <- apply(dat, 1, function(x){res <- anova(aov(x ~ conc))$`Pr(>F)`[1]})
  aov.df <- data.frame(probe = names(pvals), pval = unname(pvals), analysisID = analysisID)
  aov.df <- aov.df[aov.df$pval < 0.05, ]
  htpp.aov.filter <- rbind(htpp.aov.filter, aov.df)
}

# save ANOVA filter results
write_parquet(httr.aov.filter, paste0(data.path, "3_filtered_features/httr_aov_filter.parquet"))
write_parquet(htpp.aov.filter, paste0(data.path, "3_filtered_features/htpp_aov_filter.parquet"))


## Read in Nomic IDs and convert to httr probe IDs
httr.probes <- read.csv(paste0(data.path, "probemaps/100739_homo_sapiens_wt_probemap.csv"))
nomic <- read.csv(paste0(data.path, "probemaps/Nomic_probes_S1500_jan2024.csv"))
nomic <- nomic[nomic$Status != "On roadmap", ]
nomic <- merge(nomic, httr.probes, by = "ENTREZ_ID", all = FALSE)
nomic.probes <- unique(nomic$PROBE_NAME)

## Read in S1500 probes
s1500 <- read.csv(paste0(data.path, "probemaps/100766_homo_sapiens_s1500_probemap.csv"))
s1500.probes <- unique(s1500$PROBE_NAME)

#### Perform curve fitting for features in any list ####

# perform DR analysis for each dataset
models = c("Exp2","Exp3","Exp4","Exp5","Poly2","Lin","Power","Hill")
ncpus = 8

httr.bmd.res <- data.frame()
for(i in c(1:length(httr.paths))){
  analysisID <- httr.nms[i]
  print(analysisID)
  
  # get all IDs
  all.ids <- c(httr.wtt.filter$probe[httr.wtt.filter$analysisID == analysisID], 
               httr.aov.filter$probe[httr.aov.filter$analysisID == analysisID],
               nomic.probes, s1500.probes) %>% unique()
  
  # read in and process httr data
  lcpm <- read_parquet(httr.paths[i]) %>% as.data.frame()
  conc <- lcpm[1,-1] %>% unlist()
  rownames(lcpm) <- lcpm$probe_id
  lcpm <- lcpm[-1,-1]
  
  add.val <- abs(min(lcpm)) + 0.01*abs(min(lcpm))
  lcpm <- lcpm + add.val
  lcpm <- lcpm[rownames(lcpm) %in% all.ids, ];
  
  # perform curve fitting
  res <- PerformCurveFitting(data = lcpm, dose = conc, ncpus = ncpus, models = models)
  res <- FilterDRFit(res, lof.pval = 0)
  print("fitting done!")
  
  # perform benchmark dose calculation
  bmds <- PerformBMDCalc(res, ncpus = ncpus, num.sds = 1, sample.mean = TRUE)
  bmds <- bmds[bmds$all.pass,]
  bmds$analysisID <- analysisID
  httr.bmd.res <- rbind(httr.bmd.res, bmds)
  
}

htpp.bmd.res <- data.frame()
for(i in c(1:length(htpp.paths))){
  analysisID <- htpp.nms[i]
  print(analysisID)
  
  # read in and process httr data
  dat <- read_parquet(htpp.paths[i]) %>% as.data.frame()
  conc <- dat[1,-1] %>% unlist()
  rownames(dat) <- dat$probe_id
  dat <- dat[-1,-1]
  
  # remove features with >90% zeros
  num.zeros <- apply(dat, 1, function(x){sum(x == 0)})
  inds.remove <- which(num.zeros > 0.9*dim(dat)[2])
  dat <- dat[-inds.remove, ]
  
  add.val <- abs(min(dat)) + 0.01*abs(min(dat))
  dat <- dat + add.val
  
  # perform curve fitting
  res <- PerformCurveFitting(data = dat, dose = conc, ncpus = ncpus, models = models)
  res <- FilterDRFit(res, lof.pval = 0)
  print("fitting done!")
  
  # perform benchmark dose calculation
  bmds <- PerformBMDCalc(res, ncpus = ncpus, num.sds = 1, sample.mean = TRUE)
  bmds <- bmds[bmds$all.pass, ]
  bmds$analysisID <- analysisID
  htpp.bmd.res <- rbind(htpp.bmd.res, bmds)
}


#### Split up BMD results and save ####

# filter results by each filter type
httr.bmd.nomic <- httr.bmd.res[httr.bmd.res$gene.id %in% nomic.probes, ]
httr.bmd.s1500 <- httr.bmd.res[httr.bmd.res$gene.id %in% s1500.probes, ]

httr.bmd.aov <- data.frame()
httr.bmd.wtt <- data.frame()
for(i in c(1:length(httr.nms))){
  probes.aov <- httr.aov.filter$probe[httr.aov.filter$analysisID == httr.nms[i]]
  probes.wtt <- httr.wtt.filter$probe[httr.wtt.filter$analysisID == httr.nms[i]]
  
  temp <- httr.bmd.res[httr.bmd.res$analysisID == httr.nms[i], ]
  
  httr.bmd.aov <- rbind(httr.bmd.aov, temp[temp$gene.id %in% probes.aov, ])
  httr.bmd.wtt <- rbind(httr.bmd.wtt, temp[temp$gene.id %in% probes.wtt, ])
}

htpp.bmd.aov <- data.frame()
htpp.bmd.wtt <- data.frame()
for(i in c(1:length(htpp.nms))){
  probes.aov <- htpp.aov.filter$probe[htpp.aov.filter$analysisID == htpp.nms[i]]
  probes.wtt <- htpp.wtt.filter$probe[htpp.wtt.filter$analysisID == htpp.nms[i]]
  
  temp <- htpp.bmd.res[htpp.bmd.res$analysisID == htpp.nms[i], ]
  
  htpp.bmd.aov <- rbind(htpp.bmd.aov, temp[temp$gene.id %in% probes.aov, ])
  htpp.bmd.wtt <- rbind(htpp.bmd.wtt, temp[temp$gene.id %in% probes.wtt, ])
}

# write out all sets of BMD results
write_parquet(httr.bmd.res, paste0(data.path, "4_bmd_results/httr_bmd_all.parquet"))
write_parquet(httr.bmd.nomic, paste0(data.path, "4_bmd_results/httr_bmd_nomic.parquet"))
write_parquet(httr.bmd.s1500, paste0(data.path, "4_bmd_results/httr_bmd_s1500.parquet"))
write_parquet(httr.bmd.aov, paste0(data.path, "4_bmd_results/httr_bmd_aov.parquet"))
write_parquet(httr.bmd.wtt, paste0(data.path, "4_bmd_results/httr_bmd_wtt.parquet"))

write_parquet(htpp.bmd.res, paste0(data.path, "4_bmd_results/htpp_bmd_all.parquet"))
write_parquet(htpp.bmd.aov, paste0(data.path, "4_bmd_results/htpp_bmd_aov.parquet"))
write_parquet(htpp.bmd.wtt, paste0(data.path, "4_bmd_results/htpp_bmd_wtt.parquet"))


# fix wtt results
httr.bmds <- read_parquet(paste0(data.path, "4_bmd_results/httr_bmd_all.parquet"))
htpp.bmds <- read_parquet(paste0(data.path, "4_bmd_results/htpp_bmd_all.parquet"))

httr.wtt <- read_parquet(paste0(data.path, "3_filtered_features/httr_wtt_filter.parquet"))
htpp.wtt <- read_parquet(paste0(data.path, "3_filtered_features/htpp_wtt_filter.parquet"))



