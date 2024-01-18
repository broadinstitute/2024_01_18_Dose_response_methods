# Example BMD analysis
# Jessica Ewald
# March 10, 2021

# The curve fitting functions assume that the data matrix has already been normalized and pre-filtered
# Curves will be fit to every row in the data

source("./fastbmdR_main.R")
source("./fastbmdR_utils.R")

library(arrow)
library(dplyr)

# read in example data object
data <- read_parquet("./data.parquet") %>% as.data.frame()
rownames(data) <- data$gene_id
data <- data[,-1]

# set up parameters (model choice and associated doses)
# options are: "Exp2","Exp3","Exp4","Exp5","Poly2","Lin","Power","Hill","Poly3","Poly4"
# we usually exclude higher order polynomials (Poly3 and Poly4)
models = c("Exp2","Exp3","Exp4","Exp5","Poly2","Lin","Power","Hill")
dose = c(rep(0,5), rep(25,5), rep(100,5), rep(200,5), rep(300,5), rep(400,5))

# set the number of threads for parallel computing (defaults to 1, can be as high as your machine supports)
ncpus = 6

# perform curve fitting
start <- Sys.time()
res <- PerformCurveFitting(data = data, dose = dose, ncpus = ncpus, models = models)
end <- Sys.time()

end-start

# filter to select best fit model for each gene
res <- FilterDRFit(res, lof.pval = 0.1)

# perform benchmark dose calculation
bmd.res <- PerformBMDCalc(res, ncpus = ncpus, num.sds = 1, sample.mean = TRUE)
bmd.pass <- bmd.res[bmd.res$all.pass,]

# Results are returned for each gene in res. See the last few columns of bmd.res 
# for whether a gene BMD passed various quality criteria.

