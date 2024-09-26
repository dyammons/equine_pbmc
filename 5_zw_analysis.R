#!/usr/bin/Rscript
## This Code is for DESeq2- differential expression analysis of RNAseq gene counts
#Input_files: "~/Desktop/GJCF_PTOA_Take2/"
#Output_files: "~/Desktop/output"
#Project: GJCF, 4 horses over 5 time-points 

#Session info
sessionInfo()

## Load libraries
library(devtools)
library(Seurat)
library(clusterProfiler)
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
library(msigdbr)
library(ggplot2)

require(RColorBrewer)
require(ComplexHeatmap)
require(circlize)
require(digest)
require(cluster)
library(rmarkdown)

#load input directory and define output
genecounts <- read.csv("../input_bulk/counts.csv", row.names = 1)
head(genecounts)
ResultFile <- "../output/"

#import metadata 
# I couldn't read in you files, so I am making it manually
# the added benefit of making it like this is that you cannot accidently
# order the metadata rows differently than the columns. DEseq2 does not check
# that they match so if there is an error in ordering it will results in 
# incorrect ananlysis.
metaData <- data.frame(sample = colnames(genecounts)) %>%
  separate(sample, remove = F, sep = "_", into = c("horse", "timepoint"))

#Create the DEseq2 object; assign horse ID in columns w/ metadata in corresponding rows
dds <- DESeqDataSetFromMatrix(
  countData = genecounts, 
  colData = metaData,
  design = ~ timepoint + horse)

#view object
dds

#Filter out genes w/ counts < 10 reads and that are expressed in less than 9 horses, ideally ~50% of sample
dds <- estimateSizeFactors(dds) ### this is not needed the `DESeq` function runs it
keep <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 9 ### I wouldn't reccomend using normalized counts for filtering

#new DESeq2 object based on filtered criteria above 
dds <- dds[keep,]

#Transform and plot the data with PCA
rld <- varianceStabilizingTransformation(dds, blind=TRUE)
p <- plotPCA(rld, intgroup = "timepoint")
ggsave(paste0(ResultFile, "/pca_Day.png"), width = 7, height = 7)

p <- plotPCA(rld, intgroup = "horse")
ggsave(paste0(ResultFile, "/pca2_Horse.png"), width = 7, height = 7)

#QC heatmap
rld <- rlog(dds, blind = FALSE)
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
p <- pheatmap::pheatmap(rld_cor)
ggsave(p, file = paste0(ResultFile, "/pheatmap.png"))


#Differential expression analysis w/ DESeq2 (for more info -> ?DESeq)
dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res)

# CAUTION: The results function is doing something different than you think here.
# If you call `results` without the "name" or "contrast" argument you will only 
# get the results from the last contrast listed in the output of `resultsNames`.

resultsNames(dds)
# > resultsNames(dds)
# [1] "Intercept"                "timepoint_Day126_vs_Day0"
# [3] "timepoint_Day14_vs_Day0"  "timepoint_Day28_vs_Day0" 
# [5] "timepoint_Day56_vs_Day0"  "horse_PC02_vs_PC01"      
# [7] "horse_PC03_vs_PC01"       "horse_PC04_vs_PC01"     

# In this case that is horse_PC04_vs_PC01 which is definely not a desirable 
# contrast to evaluate how GEX is changing over the course of the study. Using 
# the Wald test (as you did here) you have to extract the results for each 
# contrast you are interested in. Below is code to extract the DEGs for the 
# contrasts that look at GEX over time. (there are other ways to get there too)

res.list <- lapply(resultsNames(dds)[2:5], function(name){
    res <- as.data.frame(results(
        dds, 
        name = name,
        lfcThreshold = 0
    )) %>%
        rownames_to_column(var = "gene") %>%
        mutate(contrast = name)
})
res <- do.call(rbind, res.list) %>%
    filter(padj < 0.05)



#reorder results based on padjust -- the code defining `sig_res` orders by padj

# The above approach works, but the rowname indexes get mixed up
# Below is a dplyr alternative

#Extract the results and save as a .csv
sig_res <- arrange(res, padj)

write.csv(
  sig_res, quote = FALSE, row.names = FALSE,
  file = paste0( ResultFile, "SDEGs.csv")
)

# I would reccomend using dplyr for all plotting (except in certain cases such
# as heatmaps)

#Sig DEG heat map
mat <- counts(dds, normalized = TRUE)[unique(sig_res$gene), ]
mat_scaled <- t(scale(t(mat)))
M <- (1 - cor(t(mat_scaled), method = "pearson")) / 2
hc <- hclust(as.dist(M), method = "complete")

#Generate heatmap
ht <- Heatmap(
  name = "Scaled\nexpression",
  mat_scaled, 
  cluster_columns = F, 
  show_row_names = F,
  row_split = factor(
    cutree(hc, k = 6),
    levels = c(4, 2, 1, 6, 5, 3)
  ),
  column_split = factor(
    gsub("PC.._", "", colnames(mat)),
    levels = unique(gsub("PC.._", "", colnames(mat)))
  )
)
png(file = paste0(ResultFile, "/DE_heat_.png"), width = 4000, height = 4000, res = 400)
par(mfcol = c(1, 1))
draw(ht)
dev.off()

#Generate boxplot of average expression at each timepoint
avg_exp <- as.data.frame(t(mat_scaled)) %>%
  rownames_to_column() %>%
  pivot_longer(cols = colnames(.)[2:ncol(.)], names_to = "gene") %>%
  separate(rowname, remove = F, sep = "_", into = c("horse", "timepoint")) %>%
  group_by(gene, timepoint) %>%
  summarize(
    MEAN = mean(value)
  )
head (avg_exp)

avg_exp$timepoint <- factor(
  avg_exp$timepoint, levels = c("Day0", "Day14", "Day28", "Day56", "Day126")
)

#Bring over the clustering results -- could improve on this
avg_exp <- left_join(
  avg_exp,
  rownames_to_column(data.frame(
    clade = cutree(hc, k = 6)
  )),
  by = c("gene" = "rowname")
)
head (avg_exp)

#Create the plot
ggplot(data = avg_exp, aes(x = timepoint, y = MEAN)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(clade ~ ., nrow = 2)
ggsave(file = paste0(ResultFile, "/gex_over_time.png"))

#Print the genes from the clade that decrease with time
print(unique(avg_exp[avg_exp$clade == 2, ]$gene))

#Save the clustering results to complete enrichment scoring w/ scrna data
write.csv(
  avg_exp, quote = FALSE, row.names = FALSE,
  file = paste0(ResultFile, "/deg_clades.csv")
)

### (3) Run GSEA of the LRT clades
avg_exp <- avg_exp[! duplicated(avg_exp$gene), ]
degs <- split(avg_exp$gene, avg_exp$clade)

#Use clusterProfiler
gene_sets <- as.data.frame(
  msigdbr(species = "horse", category = "C2", subcategory = "REACTOME")
) %>% 
  distinct(gs_name, gene_symbol)
res <- lapply(degs, enricher, TERM2GENE = gene_sets)
res <- lapply(res, as.data.frame)
res <- bind_rows(res, .id = 'clade')
res$ID <- gsub("REACTOME_", "", res$ID)

#Make another heatmap
ht <- Heatmap(
  name = "Scaled\nexpression",
  mat_scaled, 
  cluster_columns = F, 
  show_row_names = F,
  #   cluster_rows = hc
  row_split = factor(
    cutree(hc, k = 6),
    levels = c(4, 2, 1, 6, 5, 3)
  ),
  column_split = factor(
    gsub("PC.._", "", colnames(mat)),
    levels = unique(gsub("PC.._", "", colnames(mat)))
  ),
  right_annotation = rowAnnotation(
    textbox = anno_textbox(
      factor(cutree(hc, k = 6), levels = c(4, 2, 1, 6, 5, 3)), 
      split(res$ID, res$clade),
      by = "anno_block",
      gp = gpar(fontsize = 10, col = "black"),
      background_gp = gpar(fill = NA, col = NA),
      max_width = unit(2.25, "in")
    )
  )
)
png(file = paste0(ResultFile, "/DE_heat.png"), width = 4500, height = 4000, res = 400)
par(mfcol = c(1, 1))
draw(ht)
dev.off()
