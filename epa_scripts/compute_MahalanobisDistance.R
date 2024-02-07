library(arrow)
library(dplyr)
library(edgeR)

source("./globalMahalanobisDistances_v4.R")

Table1 = read_parquet("../2024_01_18_Dose_response_methods-data/data/1_orig_data/htpp_well_norm.parquet")

# try using simplified functions
meta <- Table1[,c(1:11)]
CPdata <- Table1[,-c(1:11)] %>% as.matrix()
CPdata <- CPdata[meta$plate_id == "TC00000483" & (meta$chem_id == "Actinomycin D" | meta$chem_id == "Dimethyl Sulfoxide"),!(colnames(CPdata) %in% c("f_1073", "f_1074"))] # these have constant or nearly constant values
meta <- meta[meta$plate_id == "TC00000483" & (meta$chem_id == "Actinomycin D" | meta$chem_id == "Dimethyl Sulfoxide"), ]

md.input <- prepMahalanobisDistances(CPdata, 0.95, meta$dose_level)

plate.inds <- meta$plate_id == "TC00000483"
meta.plate <- meta[plate.inds, ]
CP.plate <- CPdata[plate.inds, ]

meta$md <- computeMahalanobisDistance(CPdata, md.input$RotationMatrix, md.input$invCov, as.character(meta$dose_level), "0")
plot(meta$dose_level, 
     meta$md,
     xlab = "Dose level", 
     ylab = "global Mahalanobis Distance",
     main = "Actinomycin D: Cell Painting (single chem)")


meta.plate$md <- computeMahalanobisDistance(CP.plate, md.input$RotationMatrix, md.input$invCov, meta.plate$chem_id, "Dimethyl Sulfoxide")
plot(meta.plate$dose_level[meta.plate$chem_id == "Actinomycin D" | meta.plate$chem_id == "Dimethyl Sulfoxide"], 
     meta.plate$md[meta.plate$chem_id == "Actinomycin D" | meta.plate$chem_id == "Dimethyl Sulfoxide"],
     xlab = "Dose level", 
     ylab = "global Mahalanobis Distance",
     main = "Actinomycin D: Cell Painting")


### try with transcriptomics data
httr <- read_parquet("../2024_01_18_Dose_response_methods-data/data/1_orig_data/httr_counts.parquet") %>% as.data.frame()
httr.meta <- read_parquet("../2024_01_18_Dose_response_methods-data/data/1_orig_data/httr_key.parquet") %>% as.data.frame()

rownames(httr) <- httr$probe_id
httr <- httr[,-1]
httr <- as.matrix(httr)

httr <- httr[,colnames(httr) %in% httr.meta$sample_id[httr.meta$qc_flag == "OK"]]
httr.meta <- httr.meta[httr.meta$qc_flag == "OK", ]
httr.meta$chem_name[is.na(httr.meta$chem_name)] <- "Dimethyl Sulfoxide"

# filter based on abundance
keep.feat <- apply(httr, 1, mean) > 5
httr <- httr[keep.feat, ]

# convert to log-cpm
nf <- calcNormFactors(httr, method="RLE")
lcpm <- cpm(httr, lib.size=colSums(httr)*nf, log = TRUE, prior.count = 1)
lcpm <- t(lcpm)

inds <- httr.meta$chem_name == "Actinomycin D" | httr.meta$chem_name == "Dimethyl Sulfoxide"
httr.meta <- httr.meta[inds, ]
lcpm <- lcpm[inds, ]

md.input <- prepMahalanobisDistances(lcpm, 0.95, as.character(httr.meta$dose_level))
httr.meta$md <- computeMahalanobisDistance(lcpm, md.input$RotationMatrix, md.input$invCov, as.character(httr.meta$dose_level), "0")
plot(httr.meta$dose_level, 
     httr.meta$md,
     xlab = "Dose level", 
     ylab = "global Mahalanobis Distance",
     main = "Actinomycin D: transcriptomics")


plot(httr.meta$dose_level[httr.meta$chem_name == "Actinomycin D" | httr.meta$chem_name == "Dimethyl Sulfoxide"], 
     httr.meta$md[httr.meta$chem_name == "Actinomycin D" | httr.meta$chem_name == "Dimethyl Sulfoxide"],
     xlab = "Dose level", 
     ylab = "global Mahalanobis Distance",
     main = "Actinomycin D: transcriptomics")

