# Develop PLS-based method
# Jessica Ewald
# February 2, 2024

# setwd("/Users/jessicaewald/NetbeansProjects/2024_01_18_Dose_response_methods/analysis")

library(arrow)
library(mdatools)
library(dplyr)
source("../fastbmdR/fastbmdR_main.R")
source("../fastbmdR/fastbmdR_utils.R")

## try with my data
df <- as.data.frame(read_parquet("/Users/jessicaewald/NetbeansProjects/2024_01_18_Dose_response_methods/2024_01_18_Dose_response_methods-data/data/2_sampled_data/httr/3_reps/Dexamethasone_3reps_10_httr.parquet"))
rownames(df) <- df$probe_id
df <- df[ ,-1]
dose <- df[rownames(df) == "Dose", ] %>% unlist()
df <- df[-1, ]
df <- t(df)

dose <- data.frame(dose = dose, sample_id = rownames(df))
dose.unique <- unique(dose$dose)
dose.rank <- rank(dose.unique)-1
dose.df <- data.frame(dose = dose.unique, rank = dose.rank)
dose <- merge(dose, dose.df, by = "dose")
dose <- dose[match(rownames(df), dose$sample_id), ]

res <- computePLS(df, dose$dose, 5)
bmds <- scoresPOD(res$scores, dose$dose)

# Need to decide filters:
# all BMD quality filters must pass
# fit to regular dose
# top component with a BMD that passes 

