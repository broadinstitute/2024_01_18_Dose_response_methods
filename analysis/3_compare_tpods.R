#### Compute differences in tPOD distributions across conditions ####

# Metrics: 
# KS p-value for tPOD dists between conditions, for each of the 11 chemicals (paired?)
# Paired difference of tPOD dist mean (1 test for all 11 chemicals) between conditions
# Paired difference of tPOD dist variance (1 test for all 11 chemicals) between conditions

setwd("/Users/jessicaewald/NetbeansProjects/2024_01_18_Dose_response_methods/analysis")
library(arrow)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)

##### Read in all tpods #####

# Read in tPOD results
tpod.uni <- read_parquet("../2024_01_18_Dose_response_methods-data/data/5_tpod_results/univar_tpods.parquet")
tpod.multi.all <- read_parquet("../2024_01_18_Dose_response_methods-data/data/5_tpod_results/multivar_tpods.parquet")
tpod.multi.s1500 <- read_parquet("../2024_01_18_Dose_response_methods-data/data/5_tpod_results/multivar_tpods_s1500.parquet")
tpod.multi.nomic <- read_parquet("../2024_01_18_Dose_response_methods-data/data/5_tpod_results/multivar_tpods_nomic.parquet")

# Combine all multivar together
# all features
info <- lapply(tpod.multi.all$AnalysisID, function(x) {unlist(strsplit(x, "_"))})
info <- unlist(info)
info <- matrix(info, ncol = 4, byrow = TRUE) %>% as.data.frame()
colnames(info) <- c("Compound", "Replicates", "Iteration", "Platform")
info$Features <- "All"
info$AnalysisID <- tpod.multi.all$AnalysisID
tpod.multi <- cbind(info, tpod.multi.all[,-1])

# s1500
info <- lapply(tpod.multi.s1500$AnalysisID, function(x) {unlist(strsplit(x, "_"))})
info <- unlist(info)
info <- matrix(info, ncol = 4, byrow = TRUE) %>% as.data.frame()
colnames(info) <- c("Compound", "Replicates", "Iteration", "Platform")
info$Features <- "S1500"
info$AnalysisID <- tpod.multi.s1500$AnalysisID
tpod.temp <- cbind(info, tpod.multi.s1500[,-1])
tpod.multi <- rbind(tpod.multi, tpod.temp)

info <- lapply(tpod.multi.nomic$AnalysisID, function(x) {unlist(strsplit(x, "_"))})
info <- unlist(info)
info <- matrix(info, ncol = 4, byrow = TRUE) %>% as.data.frame()
colnames(info) <- c("Compound", "Replicates", "Iteration", "Platform")
info$Features <- "Nomic"
info$AnalysisID <- tpod.multi.nomic$AnalysisID
tpod.temp <- cbind(info, tpod.multi.nomic[,-1])
tpod.multi <- rbind(tpod.multi, tpod.temp)


# try combine everything
tpod.uni.melt <- reshape2::melt(tpod.uni, id.vars = colnames(tpod.uni)[1:6], measure.vars = colnames(tpod.uni)[7:11], variable.name = "tpod_type")
tpod.multi.melt <- reshape2::melt(tpod.multi, id.vars = colnames(tpod.multi)[1:6], measure.vars = colnames(tpod.multi)[c(7,10,13)], variable.name = "tpod_type")
tpod.melt <- rbind(tpod.uni.melt, tpod.multi.melt)

# Use log10 values
tpod.melt$log10_pod <- log10(tpod.melt$value)

# Give more identifiers
tpod.melt$unique_data <- paste0(tpod.melt$Platform, "_", tpod.melt$Features)

tpod.melt$pod_category <- NA
tpod.melt$pod_category[grepl("multi", tpod.melt$tpod_type)] <- "Multivariate"
tpod.melt$pod_category[grepl("uni", tpod.melt$tpod_type)] <- "Univariate"

tpod.melt$pod_id <- paste0("pod", c(1:dim(tpod.melt)[1]))

scenarios <- tpod.melt[,c("Compound", "Replicates", "Platform", "Features", "tpod_type")] %>% distinct()
scenarios$scenario <- paste0("scenario", c(1:dim(scenarios)[1]))
tpod.melt <- merge(tpod.melt, scenarios, by = c("Compound", "Replicates", "Platform", "Features", "tpod_type"))
tpod.melt <- tpod.melt[tpod.melt$tpod_type != "POD_multi_lmd", ]


### Now compare POD types
# restrict to 3 reps, htpp (all features) & httr (wtt filter)

pod.types <- as.data.table(tpod.melt)
pod.types <- pod.types[Replicates == "3reps" & ((Platform == "htpp" & Features == "All") | (Platform == "httr" & (Features == "WTT" | Features == "All")))]

# Rename types
pod.types$tpod_type <- as.character(pod.types$tpod_type)
pod.types$tpod_type[pod.types$tpod_type == "POD_uni_centile"] <- "10th perc."
pod.types$tpod_type[pod.types$tpod_type == "POD_uni_mode"] <- "First mode"
pod.types$tpod_type[pod.types$tpod_type == "POD_uni_lcrd"] <- "LCRD"
pod.types$tpod_type[pod.types$tpod_type == "POD_uni_gs"] <- "Gene set"
pod.types$tpod_type[pod.types$tpod_type == "POD_uni_rs"] <- "Random set"
pod.types$tpod_type[pod.types$tpod_type == "POD_multi_gmd"] <- "global MD"
pod.types$tpod_type[pod.types$tpod_type == "POD_multi_pls"] <- "PLS-DA"
pod.types$tpod_type <- factor(pod.types$tpod_type, 
                              levels = c("10th perc.", "First mode", "LCRD", "Gene set", "Random set",
                                         "global MD", "PLS-DA"))

plot.pods <- pod.types
plot.pods <- pod.types[pod.types$Compound %in% c("Actinomycin D", "Cucurbitacin I", "Cycloheximide", "Trichostatin A"), ]

new.cols <- c("#E69F00", "#56B4E9")

# make plot
bp <- ggplot(plot.pods, aes(x = log10_pod, y = Platform, color = Platform, fill = Platform)) +
  geom_boxplot(alpha=0.3) +
  scale_color_manual(values = new.cols) +
  scale_fill_manual(values = new.cols) +
  facet_grid(cols = vars(Compound), rows = vars(tpod_type), scale = "free") +
  theme_test() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        #legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

png("../2024_01_18_Dose_response_methods-data/data/poster_figures/pod_types.png",
    units = "in", res = 300, height = 7, width = 6)
print(bp)
dev.off()

ggplot(plot.pods, aes(x = log10_pod, y = Platform, fill = Platform)) +
  geom_violin() +
  scale_fill_manual(values = new.cols) +
  facet_grid(cols = vars(Compound), rows = vars(tpod_type), scale = "free") +
  theme_test() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

# compute stats for types
type.stats <- as.data.table(pod.types)
type.stats[, signif(sum(is.na(log10_pod))/.N, 2), by = .(tpod_type)]
type.stats <- type.stats[,.(pod_sd = sd(value, na.rm = T), log_sd = sd(log10_pod, na.rm = T), 
                            pod_var = var(value, na.rm = T), log_var = var(log10_pod, na.rm = T),
                            pod_mean = mean(value, na.rm = T), log_mean = mean(log10_pod, na.rm = T),
                            pod_iqr = IQR(value, na.rm = T), log_iqr = IQR(log10_pod, na.rm = T),
                            fr = signif(sum(is.na(log10_pod))/.N, 2)), 
                         by = .(tpod_type, Compound, Platform)]


# Visualize the relationship
plot(type.stats$pod_mean, type.stats$pod_sd)
plot(type.stats$pod_mean, type.stats$pod_var)
plot(type.stats$pod_mean, type.stats$pod_iqr)

plot(type.stats$log_mean, type.stats$log_sd)
plot(type.stats$log_mean, type.stats$log_var)
plot(type.stats$log_mean, type.stats$log_iqr)

stat.summary <- type.stats[, .(mean_sd = signif(mean(log_sd, na.rm = T),2), mean_var = mean(log_var, na.rm = T), mean_iqr = mean(log_iqr, na.rm = T)), 
                           by = .(tpod_type)]

#### Compute summary stats ####

# compute summary statistics
tpod.stats <- as.data.table(tpod.melt)
tpod.stats <- tpod.stats[,.(Compound = Compound, Replicates = Replicates, Platform = Platform, Features = Features, POD_type = tpod_type, 
                            POD_mean = mean(log10_pod, na.rm = TRUE), 
                            POD_median = median(log10_pod, na.rm = TRUE),
                            POD_sd = var(log10_pod, na.rm = TRUE),
                            POD_mad = mad(log10_pod, na.rm = TRUE),
                            POD_iqr = IQR(log10_pod, na.rm = TRUE)), 
                         by = .(scenario, Compound)] %>% distinct() %>% na.omit()

# Compare 2 reps vs 3 reps
# Conclusion: mean POD unchanged, but variance slightly decreased for 3 replicates
tpod.reps <- tpod.stats
tpod.reps$id_col <- paste0(tpod.reps$Compound, "_", tpod.reps$Platform, "_", tpod.reps$Features, "_", tpod.reps$POD_type)

tpod.reps.mean <- dcast(tpod.reps, id_col ~ Replicates, value.var = "POD_mean")
tpod.reps.mean$diff <- tpod.reps.mean[,"2reps"] - tpod.reps.mean[, "3reps"]
tpod.reps.mean$fc <- (10^tpod.reps.mean[,"2reps"])/(10^tpod.reps.mean[, "3reps"])
mean(tpod.reps.mean$fc)
t.test(x = unlist(tpod.reps.mean[,"2reps"]), y = unlist(tpod.reps.mean[,"3reps"]), paired = T)
# no significant difference

tpod.reps.iqr <- dcast(tpod.reps, id_col ~ Replicates, value.var = "POD_iqr")
tpod.reps.iqr$diff <- tpod.reps.iqr[,"2reps"] - tpod.reps.iqr[,"3reps"] 
tpod.reps.iqr$fc <- (10^tpod.reps.iqr[,"2reps"])/(10^tpod.reps.iqr[, "3reps"])
t.test(x = unlist(tpod.reps.iqr[,"2reps"]), y = unlist(tpod.reps.iqr[,"3reps"]), paired = T)
mean(tpod.reps.iqr$fc)
# p = 7.4 e-12


# Compare platform
# Conclusion: HTPP has on average a lower POD. Variance is centered around zero. 
tpod.platform <- tpod.stats
tpod.platform$id_col <- paste0(tpod.platform$Compound, "_", tpod.platform$Replicates, "_", tpod.platform$Features, "_", tpod.platform$POD_type)

tpod.platform.mean <- dcast(tpod.platform, id_col ~ Platform, value.var = "POD_mean") %>% na.omit()
tpod.platform.mean$diff <- tpod.platform.mean[,"htpp"] - tpod.platform.mean[,"httr"]
tpod.platform.mean$fc <- (10^tpod.platform.mean$htpp)/(10^tpod.platform.mean$httr)
mean(tpod.platform.mean$fc)
t.test(x = unlist(tpod.platform.mean[,"htpp"]), y = unlist(tpod.platform.mean[,"httr"]), paired = T)
# htpp ~10-fold lower, p = 0.0011

tpod.platform.iqr <- dcast(tpod.platform, id_col ~ Platform, value.var = "POD_iqr") %>% na.omit()
tpod.platform.iqr$diff <- tpod.platform.iqr[,"htpp"] - tpod.platform.iqr[,"httr"]
t.test(x = unlist(tpod.platform.iqr[,"htpp"]), y = unlist(tpod.platform.iqr[,"httr"]), paired = T)
mean(tpod.platform.iqr$diff)
# SD of Cell Painting was lower than transcriptomics (mean difference = 0.87 SDs). P-val = 0.0044.

# Compare WTT vs. ANOVA
tpod.filter <- tpod.stats[Features %in% c("WTT", "ANOVA")]
tpod.filter$id_col <- paste0(tpod.filter$Compound, "_", tpod.filter$Replicates, "_", tpod.filter$Platform, "_", tpod.filter$POD_type)

tpod.filter.mean <- dcast(tpod.filter, id_col ~ Features, value.var = "POD_mean") %>% na.omit()
tpod.filter.mean$diff <- tpod.filter.mean[,"ANOVA"] - tpod.filter.mean[,"WTT"]
tpod.filter.mean$fc <- (10^tpod.filter.mean$ANOVA)/(10^tpod.filter.mean$WTT)
t.test(x = unlist(tpod.filter.mean[,"ANOVA"]), y = unlist(tpod.filter.mean[,"WTT"]), paired = T)
# ns diff

tpod.filter.iqr <- dcast(tpod.filter, id_col ~ Features, value.var = "POD_iqr") %>% na.omit()
tpod.filter.iqr$diff <- tpod.filter.iqr[,"ANOVA"] - tpod.filter.iqr[,"WTT"]
t.test(x = unlist(tpod.filter.iqr[,"ANOVA"]), y = unlist(tpod.filter.iqr[,"WTT"]), paired = T)
# ns diff

# Compare S1500 vs. All (WTT) - httr only
tpod.s1500 <- tpod.stats[Features %in% c("WTT", "S1500") & Platform == "httr"]
tpod.s1500$id_col <- paste0(tpod.s1500$Compound, "_", tpod.s1500$Replicates, "_", tpod.s1500$POD_type)

tpod.s1500.mean <- dcast(tpod.s1500, id_col ~ Features, value.var = "POD_mean") %>% na.omit()
tpod.s1500.mean$fc <- (10^tpod.s1500.mean$S1500)/(10^tpod.s1500.mean$WTT)
t.test(x = unlist(tpod.s1500.mean[,"S1500"]), y = unlist(tpod.s1500.mean[,"WTT"]), paired = T)


tpod.s1500.iqr <- dcast(tpod.s1500, id_col ~ Features, value.var = "POD_iqr") %>% na.omit()
tpod.s1500.iqr$diff <- tpod.s1500.iqr[,2] - tpod.s1500.iqr[,3]
tpod.s1500.iqr$fc <- (10^tpod.s1500.iqr[,"S1500"])/(10^tpod.s1500.iqr[, "WTT"])
mean(tpod.s1500.iqr$fc)
t.test(x = unlist(tpod.s1500.iqr[,"S1500"]), y = unlist(tpod.s1500.iqr[,"WTT"]), paired = T)


# Compare Nomic vs. All (WTT) - httr only
tpod.nomic <- tpod.stats[Features %in% c("WTT", "Nomic") & Platform == "httr"]
tpod.nomic$id_col <- paste0(tpod.nomic$Compound, "_", tpod.nomic$Replicates, "_", tpod.nomic$POD_type)

tpod.nomic.mean <- dcast(tpod.nomic, id_col ~ Features, value.var = "POD_mean") %>% na.omit()
tpod.nomic.mean$diff <- tpod.nomic.mean[,2] - tpod.nomic.mean[,3]
hist(tpod.nomic.mean$diff)

tpod.nomic.sd <- dcast(tpod.nomic, id_col ~ Features, value.var = "POD_sd") %>% na.omit()
tpod.nomic.sd$diff <- tpod.nomic.sd[,2] - tpod.nomic.sd[,3]
hist(tpod.nomic.sd$diff)


## POD density for poster
ggplot(tpod.melt[tpod.melt$scenario == "scenario170", ], aes(x = log10_pod)) +
  geom_histogram(bins=10) +
  theme_test() +
  ylab("Count") +
  xlab("POD") +
  theme(text=element_text(size=21))


### Analyze FC tpods
fc_ids = c("FC1", "FC2", "FC3", "FC4", "FC5")
fc_tpods = read_parquet("../2024_01_18_Dose_response_methods-data/data/5_tpod_results/univar_tpods.parquet") %>% as.data.frame()
fc_tpods = fc_tpods[fc_tpods$Features %in% fc_ids, ]

fc_tpods = reshape2::melt(fc_tpods, id.vars = colnames(tpod.uni)[1:6], measure.vars = colnames(tpod.uni)[7:11], variable.name = "tpod_type")
fc_tpods$log10_pod <- log10(fc_tpods$value)

ggplot(fc_tpods[
  (fc_tpods$Replicates == "2reps") &
    (fc_tpods$tpod_type %in% c("POD_uni_centile", "POD_uni_mode", "POD_uni_lcrd", "POD_uni_gs")), ], aes(x = log10_pod, y = Features, color = Features, fill = Features)) +
  geom_boxplot(alpha=0.3) +
  facet_grid(cols = vars(Compound), rows = vars(tpod_type), scale = "free") +
  theme_test() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# summary stats
scenarios <- fc_tpods[,c("Compound", "Replicates", "Features", "tpod_type")] %>% distinct()
scenarios$scenario <- paste0("scenario", c(1:dim(scenarios)[1]))
fc_tpods <- merge(fc_tpods, scenarios, by = c("Compound", "Replicates", "Features", "tpod_type"))

tpod.stats <- as.data.table(fc_tpods)
tpod.stats <- tpod.stats[,.(Compound = Compound, Replicates = Replicates, Features = Features, POD_type = tpod_type, 
                            POD_mean = mean(log10_pod, na.rm = TRUE), 
                            POD_median = median(log10_pod, na.rm = TRUE),
                            POD_sd = var(log10_pod, na.rm = TRUE),
                            POD_mad = mad(log10_pod, na.rm = TRUE),
                            POD_iqr = IQR(log10_pod, na.rm = TRUE)), 
                         by = .(scenario, Compound)] %>% distinct() %>% na.omit()

# FC2 vs. FC1
tpod.fc <- tpod.stats[tpod.stats$Features %in% c("FC1", "FC2"), ]
tpod.fc$id_col <- paste0(tpod.fc$Compound, "_", tpod.fc$Replicates, "_", tpod.fc$POD_type)

tpod.fc.mean <- dcast(tpod.fc, id_col ~ Features, value.var = "POD_mean") %>% na.omit()
tpod.fc.mean$diff <- tpod.fc.mean[,"FC2"] - tpod.fc.mean[,"FC1"]
tpod.fc.mean$fc <- (10^tpod.fc.mean$FC2)/(10^tpod.fc.mean$FC1)
mean(tpod.fc.mean$fc)
t.test(x = unlist(tpod.fc.mean[,"FC1"]), y = unlist(tpod.fc.mean[,"FC2"]), paired = T)
# FC2 ~1.5 fold higher, p = 0.1049

tpod.fc.iqr <- dcast(tpod.fc, id_col ~ Features, value.var = "POD_iqr") %>% na.omit()
tpod.fc.iqr$diff <- tpod.fc.iqr[,"FC2"] - tpod.fc.iqr[,"FC1"]
t.test(x = unlist(tpod.fc.iqr[,"FC2"]), y = unlist(tpod.fc.iqr[,"FC1"]), paired = T)
mean(tpod.fc.iqr$diff)
# IQR 0.05 wider, p = 0.13


# FC3 vs. FC2
tpod.fc <- tpod.stats[tpod.stats$Features %in% c("FC2", "FC3"), ]
tpod.fc$id_col <- paste0(tpod.fc$Compound, "_", tpod.fc$Replicates, "_", tpod.fc$POD_type)

tpod.fc.mean <- dcast(tpod.fc, id_col ~ Features, value.var = "POD_mean") %>% na.omit()
tpod.fc.mean$diff <- tpod.fc.mean[,"FC3"] - tpod.fc.mean[,"FC2"]
tpod.fc.mean$fc <- (10^tpod.fc.mean$FC3)/(10^tpod.fc.mean$FC2)
mean(tpod.fc.mean$fc)
t.test(x = unlist(tpod.fc.mean[,"FC2"]), y = unlist(tpod.fc.mean[,"FC3"]), paired = T)
# FC3 ~3 fold higher, p = 1.004e-05

tpod.fc.iqr <- dcast(tpod.fc, id_col ~ Features, value.var = "POD_iqr") %>% na.omit()
tpod.fc.iqr$diff <- tpod.fc.iqr[,"FC3"] - tpod.fc.iqr[,"FC2"]
t.test(x = unlist(tpod.fc.iqr[,"FC3"]), y = unlist(tpod.fc.iqr[,"FC2"]), paired = T)
mean(tpod.fc.iqr$diff)
# IQR 0.03 wider, p = 0.30


# FC4 vs. FC3
tpod.fc <- tpod.stats[tpod.stats$Features %in% c("FC3", "FC4"), ]
tpod.fc$id_col <- paste0(tpod.fc$Compound, "_", tpod.fc$Replicates, "_", tpod.fc$POD_type)

tpod.fc.mean <- dcast(tpod.fc, id_col ~ Features, value.var = "POD_mean") %>% na.omit()
tpod.fc.mean$diff <- tpod.fc.mean[,"FC4"] - tpod.fc.mean[,"FC3"]
tpod.fc.mean$fc <- (10^tpod.fc.mean$FC4)/(10^tpod.fc.mean$FC3)
mean(tpod.fc.mean$fc)
t.test(x = unlist(tpod.fc.mean[,"FC3"]), y = unlist(tpod.fc.mean[,"FC4"]), paired = T)
# FC4 ~1.5 fold higher, p = 1.437e-05

tpod.fc.iqr <- dcast(tpod.fc, id_col ~ Features, value.var = "POD_iqr") %>% na.omit()
tpod.fc.iqr$diff <- tpod.fc.iqr[,"FC4"] - tpod.fc.iqr[,"FC3"]
t.test(x = unlist(tpod.fc.iqr[,"FC4"]), y = unlist(tpod.fc.iqr[,"FC3"]), paired = T)
mean(tpod.fc.iqr$diff)
# IQR 0.01 wider, p = 0.76


# FC5 vs. FC4
tpod.fc <- tpod.stats[tpod.stats$Features %in% c("FC4", "FC5"), ]
tpod.fc$id_col <- paste0(tpod.fc$Compound, "_", tpod.fc$Replicates, "_", tpod.fc$POD_type)

tpod.fc.mean <- dcast(tpod.fc, id_col ~ Features, value.var = "POD_mean") %>% na.omit()
tpod.fc.mean$diff <- tpod.fc.mean[,"FC5"] - tpod.fc.mean[,"FC4"]
tpod.fc.mean$fc <- (10^tpod.fc.mean$FC5)/(10^tpod.fc.mean$FC4)
mean(tpod.fc.mean$fc)
t.test(x = unlist(tpod.fc.mean[,"FC4"]), y = unlist(tpod.fc.mean[,"FC5"]), paired = T)
# FC4 ~1.5 fold higher, p = 0.0001214

tpod.fc.iqr <- dcast(tpod.fc, id_col ~ Features, value.var = "POD_iqr") %>% na.omit()
tpod.fc.iqr$diff <- tpod.fc.iqr[,"FC5"] - tpod.fc.iqr[,"FC4"]
t.test(x = unlist(tpod.fc.iqr[,"FC5"]), y = unlist(tpod.fc.iqr[,"FC4"]), paired = T)
mean(tpod.fc.iqr$diff)
# IQR 0.01 wider, p = 0.76