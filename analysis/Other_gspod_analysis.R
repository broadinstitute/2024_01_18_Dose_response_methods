# Compare real gene sets to random gene sets
# Jessica Ewald
# January 16, 2024



# Previous work has shown that different pathway libaries return different gene
# set PODs. I aim to see if 1) these differences are significant across 
# subsampling, and if so, 2) if the differences are explained mainly by library 
# size. 
#
# I will do this by comparing PODs across GO BP, Reactome, and KEGG and by 
# generating random libaries with the same size characteristics of each of these,
# computing PODs for all library sets, and comparing.



library(dplyr)
source("/Users/jessicaewald/Desktop/OASIS/analysis/POD_subsampling/rscripts/fastbmdR/fastbmdR_main.R")
source("/Users/jessicaewald/Desktop/OASIS/analysis/POD_subsampling/rscripts/fastbmdR/fastbmdR_utils.R")

# Read in BMD results (ANOVA, 3 relicates)
bmd.res <- readRDS("/Users/jessicaewald/Desktop/OASIS/analysis/POD_subsampling/analysis_results/1_prefilter_compare/aovBmdRes.rds")
bmd.res <- bmd.res[grep("3reps", names(bmd.res))]

all.conc <- readRDS("/Users/jessicaewald/Desktop/OASIS/analysis/POD_subsampling/analysis_results/1_prefilter_compare/allConc.rds")
all.conc <- all.conc[grep("3reps", names(all.conc))]

# Process existing pathway libraries
wt.probes <- read.csv("/Users/jessicaewald/Desktop/OASIS/analysis/POD_subsampling/data/100739_homo_sapiens_wt_probemap.csv")
universe <- unique(wt.probes$ENTREZ_ID) %>% as.character()
universe <- na.omit(universe)

gobp <- readRDS("/Users/jessicaewald/Desktop/OASIS/analysis/POD_subsampling/data/gs_libraries/go_bp.rds")
kegg <- readRDS("/Users/jessicaewald/Desktop/OASIS/analysis/POD_subsampling/data/gs_libraries/kegg.rds")
reactome <- readRDS("/Users/jessicaewald/Desktop/OASIS/analysis/POD_subsampling/data/gs_libraries/reactome.rds")

gobp.names <- data.frame(name = gobp$term, id = names(gobp$sets))
kegg.names <- data.frame(name = kegg$term, id = names(kegg$sets))
reactome.names <- data.frame(name = reactome$term, id = names(reactome$sets))

gobp.lib <- gobp[["sets"]]
kegg.lib <- kegg[["sets"]]
reactome.lib <- reactome[["sets"]]

# Create random libraries
gobp.lib.random <- lapply(gobp.lib, function(x){
  sample(universe, length(x), replace = FALSE)
})

kegg.lib.random <- lapply(kegg.lib, function(x){
  sample(universe, length(x), replace = FALSE)
})

reactome.lib.random <- lapply(reactome.lib, function(x){
  sample(universe, length(x), replace = FALSE)
})

# Define pathway POD function
pwPOD <- function(gs.lib, universe, bmds) {
  
  # hypergeometic test
  results.test <- multiHyperGeoTest(gs.lib, universe, bmds)
  results <- as.data.frame(results.test$results)
  results$id <- rownames(results)
  rownames(results) <- NULL
  
  # get BMDs from each enriched gene set
  gs.id <- as.list(results$id)
  pw.bmds <- lapply(gs.id, function(x){bmd.temp[bmd.temp$entrez %in% gs.lib[[x]], "bmd"]})
  
  # median BMD for each gene set
  bmd.med <- lapply(pw.bmds, median)
  results$bmd.med <- unlist(bmd.med)
  
  # min 3 overlap and 5%
  results$percent.overlap <- results$`Observed Hits`/results$`Gene Set Size`
  results <- results[results$`Observed Hits` > 2 & results$percent.overlap > 0.05, ]
  
  # pathway POD is lowest median BMD that passes all criteria
  return(min(results$bmd.med))
}

# Compute PODs
exp.names <- names(bmd.res)

gobp.pod <- list()
kegg.pod <- list()
reactome.pod <- list()
gobp.random.pod <- list()
kegg.random.pod <- list()
reactome.random.pod <- list()

for(i in c(1:length(exp.names))){
  print(i)
  
  chem <- gsub("_.*", "", exp.names[i])
  doses <- all.conc[[exp.names[i]]]
  
  # convert probe IDs to Entrez
  bmd.temp <- bmd.res[[exp.names[i]]]
  bmd.temp$probe.id <- gsub(".*_", "", bmd.temp$gene.id)
  inds <- match(bmd.temp$probe.id, wt.probes$PROBE_ID)
  bmd.temp$entrez <- wt.probes$ENTREZ_ID[inds] %>% as.character()
  
  # compute PODs
  gobp.pod[[i]] <- pwPOD(gobp.lib, universe, na.omit(bmd.temp$entrez))
  kegg.pod[[i]] <- pwPOD(kegg.lib, universe, na.omit(bmd.temp$entrez))
  reactome.pod[[i]] <- pwPOD(reactome.lib, universe, na.omit(bmd.temp$entrez))
  gobp.random.pod[[i]] <- pwPOD(gobp.lib.random, universe, na.omit(bmd.temp$entrez))
  kegg.random.pod[[i]] <- pwPOD(kegg.lib.random, universe, na.omit(bmd.temp$entrez))
  reactome.random.pod[[i]] <- pwPOD(reactome.lib.random, universe, na.omit(bmd.temp$entrez))
  
}
names(gobp.pod) <- exp.names
names(kegg.pod) <- exp.names
names(reactome.pod) <- exp.names
names(gobp.random.pod) <- exp.names
names(kegg.random.pod) <- exp.names
names(reactome.random.pod) <- exp.names


# Combine PODs into dataframe
pod1 <- reshape2::melt(gobp.pod)[,c(2,1)]
pod2 <- reshape2::melt(kegg.pod)[,c(2,1)]
pod3 <- reshape2::melt(reactome.pod)[,c(2,1)]
pod4 <- reshape2::melt(gobp.random.pod)[,c(2,1)]
pod5 <- reshape2::melt(kegg.random.pod)[,c(2,1)]
pod6 <- reshape2::melt(reactome.random.pod)[,c(2,1)]

pod1$type <- "gobp"
pod2$type <- "kegg"
pod3$type <- "reactome"
pod4$type <- "gobp.random"
pod5$type <- "kegg.random"
pod6$type <- "reactome.random"

pods <- rbind(pod1, pod2)
pods <- rbind(pods, pod3)
pods <- rbind(pods, pod4)
pods <- rbind(pods, pod5)
pods <- rbind(pods, pod6)

colnames(pods) <- c("dataset", "pod", "type")
pods$chem <- gsub("_.*", "", pods$dataset)
pods$pod[is.infinite(pods$pod)] <- NA
saveRDS(pods, "/Users/jessicaewald/Desktop/OASIS/analysis/POD_subsampling/analysis_results/1_prefilter_compare/pods_gs_011624.rds")

# Visualize results
library(ggplot2)
library(ggpubr)

# make each plot with individual xlim
chem.conc <- all.conc[seq(1:11)*30]
names(chem.conc) <- gsub("_.*", "", names(chem.conc))

# across three regular gene set libraries
plots <- list()
for(i in c(1:length(chem.conc))){
  plots[[i]] <- ggplot(pods[pods$type %in% c("kegg", "reactome", "gobp") & pods$chem == names(chem.conc)[i], ], 
                       aes(x = pod, fill = type, color = type)) +
    geom_vline(xintercept = chem.conc[[i]], color = "grey", alpha = 0.5) +
    geom_density(alpha=0.3) +
    scale_x_continuous(breaks = chem.conc[[i]], trans = "log10", limits = c(chem.conc[[i]][2]/10, max(chem.conc[[i]])))+
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    ggtitle(names(chem.conc)[i]) +
    theme_classic() +
    theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
names(plots) <- names(chem.conc)
plot <- ggarrange(plotlist = plots, ncol = 4, nrow = 3, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))
annotate_figure(plot, top = text_grob("Compare pathway POD between libraries", face = "bold", size = 18))

# each regular library vs its random version

# GO BP
plots <- list()
for(i in c(1:length(chem.conc))){
  plots[[i]] <- ggplot(pods[pods$type %in% c("gobp", "gobp.random") & pods$chem == names(chem.conc)[i], ], 
                       aes(x = pod, fill = type, color = type)) +
    geom_vline(xintercept = chem.conc[[i]], color = "grey", alpha = 0.5) +
    geom_density(alpha=0.3) +
    scale_x_continuous(breaks = chem.conc[[i]], trans = "log10", limits = c(chem.conc[[i]][2]/10, max(chem.conc[[i]])))+
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    ggtitle(names(chem.conc)[i]) +
    theme_classic() +
    theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
names(plots) <- names(chem.conc)
plot <- ggarrange(plotlist = plots, ncol = 4, nrow = 3, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))
annotate_figure(plot, top = text_grob("Compare GO BP to random sets", face = "bold", size = 18))


# KEGG
plots <- list()
for(i in c(1:length(chem.conc))){
  plots[[i]] <- ggplot(pods[pods$type %in% c("kegg", "kegg.random") & pods$chem == names(chem.conc)[i], ], 
                       aes(x = pod, fill = type, color = type)) +
    geom_vline(xintercept = chem.conc[[i]], color = "grey", alpha = 0.5) +
    geom_density(alpha=0.3) +
    scale_x_continuous(breaks = chem.conc[[i]], trans = "log10", limits = c(chem.conc[[i]][2]/10, max(chem.conc[[i]])))+
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    ggtitle(names(chem.conc)[i]) +
    theme_classic() +
    theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
names(plots) <- names(chem.conc)
plot <- ggarrange(plotlist = plots, ncol = 4, nrow = 3, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))
annotate_figure(plot, top = text_grob("Compare KEGG to random sets", face = "bold", size = 18))


# Reactome
plots <- list()
for(i in c(1:length(chem.conc))){
  plots[[i]] <- ggplot(pods[pods$type %in% c("reactome", "reactome.random") & pods$chem == names(chem.conc)[i], ], 
                       aes(x = pod, fill = type, color = type)) +
    geom_vline(xintercept = chem.conc[[i]], color = "grey", alpha = 0.5) +
    geom_density(alpha=0.3) +
    scale_x_continuous(breaks = chem.conc[[i]], trans = "log10", limits = c(chem.conc[[i]][2]/10, max(chem.conc[[i]])))+
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    ggtitle(names(chem.conc)[i]) +
    theme_classic() +
    theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
names(plots) <- names(chem.conc)
plot <- ggarrange(plotlist = plots, ncol = 4, nrow = 3, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))
annotate_figure(plot, top = text_grob("Compare Reactome to random sets", face = "bold", size = 18))


# compute statistics between POD distributions
pods <- readRDS("/Users/jessicaewald/Desktop/OASIS/analysis/POD_subsampling/analysis_results/1_prefilter_compare/pods_gs_011624.rds")
chems <- names(chem.conc)


# GO BP
for(i in c(1:length(chems))){
  res <- ks.test(x = pods$pod[pods$type == "gobp" & pods$chem == chems[i]],
                 y = pods$pod[pods$type == "gobp.random" & pods$chem == chems[i]],
                 alternative = "two.sided")
  print(paste0(chems[i], ": ", res$p.value))
}

# KEGG
for(i in c(1:length(chems))){
  res <- ks.test(x = pods$pod[pods$type == "kegg" & pods$chem == chems[i]],
                 y = pods$pod[pods$type == "kegg.random" & pods$chem == chems[i]],
                 alternative = "two.sided")
  print(paste0(chems[i], ": ", res$p.value))
}

# Reactome
for(i in c(1:length(chems))){
  res <- ks.test(x = pods$pod[pods$type == "reactome" & pods$chem == chems[i]],
                 y = pods$pod[pods$type == "reactome.random" & pods$chem == chems[i]],
                 alternative = "two.sided")
  print(paste0(chems[i], ": ", res$p.value))
}



