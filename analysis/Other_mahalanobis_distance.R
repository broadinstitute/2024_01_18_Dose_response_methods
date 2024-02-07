# Develop mahalanobis metric
# Jessica Ewald
# Feb 2, 2024

#### Compute Mahalanobis distances for dose matrix
library(arrow)
library(dplyr)
library(ggplot2)
library(edgeR)

df <- as.data.frame(read_parquet("/Users/jessicaewald/NetbeansProjects/2024_01_18_Dose_response_methods/2024_01_18_Dose_response_methods-data/data/2_sampled_data/httr/3_reps/Actinomycin D_3reps_1.parquet"))
rownames(df) <- df$probe_id
df <- df[ ,-1]
dose <- df[rownames(df) == "Dose", ] %>% unlist()
df <- df[-1, ]

mahalDist <- function(df, dose){
  
  df <- t(df)
  
  # project all samples of PCs of controls
  pr.control <- prcomp(df[dose == 0, ], center = TRUE, scale = TRUE)
  df.proj <- scale(df, pr.control$center, pr.control$scale) %*% pr.control$rotation
  
  # compute Mahalanobis distance to center of controls
  df.center <- colMeans(pr.control$x)
  df.cov <- cov(pr.control$x)
  df.dist <- mahalanobis(df.proj, center = df.center, cov = df.cov, tol=1e-50)
  
}


# Conclusion: from PCA, looks like there are some pretty major batch effects? 
# The mahalanobis distances are crazy and don't make any sense.
x0 <- df.proj %>% as.data.frame()
x0$dose = dose
x0$dose[x0$dose == 0] <- sort(unique(dose))[2]/(sort(unique(dose))[3]/sort(unique(dose))[2])
x0$dose <- log10(x0$dose)

ggplot(x0, aes(x = PC1, y = PC2, color = dose)) +
  geom_point() +
  theme_bw()


# Check with full httr matrix
httr <- read_parquet("./2024_01_18_Dose_response_methods-data/data/1_orig_data/httr_counts.parquet") %>% as.data.frame()
meta <- read_parquet("./2024_01_18_Dose_response_methods-data/data/1_orig_data/httr_key.parquet") %>% as.data.frame()

rownames(httr) <- httr$probe_id
httr <- httr[,-1]
httr <- as.matrix(httr)

# add register to meta
reg.let.1 <- c("A", "C", "E", "G", "I", "K", "M", "O")
reg.let.2 <- c("B", "D", "F", "H", "J", "L", "N", "P")

reg.num.1 <- c("01", "03", "05", "07", "09", "11", "13", "15", "17", "19", "21", "23")
reg.num.2 <- c("02", "04", "06", "08", "10", "12", "14", "16", "18", "20", "22", "24")

reg1 <- lapply(reg.let.1, function(x){paste0(x, reg.num.1)}) %>% unlist()
reg2 <- lapply(reg.let.1, function(x){paste0(x, reg.num.2)}) %>% unlist()
reg3 <- lapply(reg.let.2, function(x){paste0(x, reg.num.1)}) %>% unlist()
reg4 <- lapply(reg.let.2, function(x){paste0(x, reg.num.2)}) %>% unlist()

meta$register = NA
meta$register[meta$well_id %in% reg1] <- "reg1"
meta$register[meta$well_id %in% reg2] <- "reg2"
meta$register[meta$well_id %in% reg3] <- "reg3"
meta$register[meta$well_id %in% reg4] <- "reg4"

meta$reg_plate <- paste0(meta$plate_id, "_", meta$register)

httr <- httr[,colnames(httr) %in% meta$sample_id[meta$qc_flag == "OK"]]
meta <- meta[meta$qc_flag == "OK", ]

# filter based on abundance
keep.feat <- apply(httr, 1, mean) > 5
httr <- httr[keep.feat, ]

# convert to log-cpm
nf <- calcNormFactors(httr, method="RLE")
lcpm <- cpm(httr, lib.size=colSums(httr)*nf, log = TRUE, prior.count = 1)

pr <- prcomp(t(lcpm), center = T, scale = T)
x <- pr$x %>% as.data.frame()
x <- merge(x, meta, by.x = "row.names", by.y = "sample_id")

ggplot(x, aes(x = PC1, y = PC2, color = stype)) +
  geom_point() +
  theme_bw() +
  ggtitle("Transcriptomics data")

ggplot(x[x$stype == "test sample", ], aes(x = PC1, y = PC2, color = chem_name)) +
  geom_point() +
  theme_bw()

ggplot(x, aes(x = PC1, y = PC2, color = register)) +
  geom_point() +
  theme_bw() +
  ggtitle("Transcriptomics data")

ggplot(x, aes(x = PC1, y = PC2, color = plate_id)) +
  geom_point() +
  theme_bw() +
  ggtitle("Transcriptomics data")

ggplot(x, aes(x = PC1, y = PC2, color = reg_plate)) +
  geom_point() +
  theme_bw() +
  ggtitle("Transcriptomics data")

# Check out htpp data
htpp <- read_parquet("/Users/jessicaewald/NetbeansProjects/2024_01_18_Dose_response_methods/2024_01_18_Dose_response_methods-data/data/1_orig_data/htpp_well_norm.parquet") %>% as.data.frame()
htpp.meta <- htpp[,c(1:11)]
htpp <- htpp[,-c(1:11)] %>% t()
colnames(htpp) <- htpp.meta$sample_id

pr2 <- prcomp(t(htpp), center = F, scale = F) # already scaled
x2 <- pr2$x %>% as.data.frame()
x2 <- merge(x2, htpp.meta, by.x = "row.names", by.y = "sample_id")

ggplot(x2, aes(x = PC1, y = PC2, color = stype)) +
  geom_point() +
  theme_bw() +
  ggtitle("Cell Painting data")

ggplot(x2[x2$stype == "viability positive control", ], aes(x = PC1, y = PC2, color = dose_level)) +
  geom_point() +
  theme_bw()

ggplot(x2, aes(x = PC1, y = PC2, color = chem_id)) +
  geom_point() +
  theme_bw()
