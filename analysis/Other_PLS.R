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
df <- as.data.frame(read_parquet("/Users/jessicaewald/NetbeansProjects/2024_01_18_Dose_response_methods/2024_01_18_Dose_response_methods-data/data/2_sampled_data/httr/3_reps/Actinomycin D_3reps_3.parquet"))
#df <- as.data.frame(read_parquet("/Users/jessicaewald/NetbeansProjects/2024_01_18_Dose_response_methods/2024_01_18_Dose_response_methods-data/data/2_sampled_data/httr/3_reps/Retinoic acid_3reps_13.parquet"))
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

res <- mdatools::pls(df, dose$rank, 10, cv = list("ven", 3), scale = T)
plot(res)

# this is for computing only a simple PLS
res <- mdatools::pls.simpls(x, y, ncomp, cv = FALSE)

# extract loadings
# xload <- res$xloadings

# extract scores
comps <- res$res$cal$xdecomp
xscores <- comps$scores
comps$expvar # amount explained in original dataset

# only include a component if more than an x% decrease compared to the previous one
# maybe 2%? 5%?
# cv.rmse <- c(res$res$cv$rmse)
# ((cv.rmse[2] - cv.rmse[1])/cv.rmse[1])*100
# ((cv.rmse[3] - cv.rmse[2])/cv.rmse[2])*100
# ((cv.rmse[4] - cv.rmse[3])/cv.rmse[3])*100
# ((cv.rmse[5] - cv.rmse[4])/cv.rmse[4])*100
# ((cv.rmse[6] - cv.rmse[5])/cv.rmse[5])*100

# try curve fitting on scores
xscores <- t(xscores)
xscores <- xscores + abs(min(xscores))*1.05 # shift above zero

models = c("Exp2","Exp3","Exp4","Exp5","Poly2","Lin","Power","Hill")

# perform curve fitting
curve.res <- PerformCurveFitting(data = xscores, dose = dose$dose, ncpus = 1, models = models)
curve.res <- FilterDRFit(curve.res, lof.pval = 0)

# perform benchmark dose calculation
bmds <- PerformBMDCalc(curve.res, ncpus = 1, num.sds = 1, sample.mean = TRUE)

plot(dose$rank, xscores[1,])
plot(dose$rank, xscores[2,])
plot(dose$rank, xscores[3,])
plot(dose$rank, xscores[4,])
plot(dose$rank, xscores[5,])



