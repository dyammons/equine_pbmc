#!/usr/bin/Rscript
## This Code is for DESeq2- differential expression analysis of RNAseq gene counts
#Input_files: "~/Desktop/GJCF_PTOA_Take2/"
#Output_files: "~/Desktop/output"
#Project: GJCF, 4 horses over 5 time-points 

#Session info
sessionInfo()

##Install Packages 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("ComplexHeatmap")
BiocManager::install("clusterProfiler")
install.packages("tidyverse")
install.packages("msigdbr")
install.packages('Seurat')
install.packages("rmarkdown")

## Load libraries
library(devtools)
install_github("jokergoo/ComplexHeatmap")
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
genecounts <- read.csv("~/Desktop/GJCF_PTOA_Take2/HTSeq_Gene_counts_copy.csv", row.names = 1)
head(genecounts)
ResultFile <- "~/Desktop/GJCF_PTOA_Take2"

#import metadata
metaData <- read.csv("~/Desktop/GJCF_PTOA_Take2/GJCF_Metadata.csv", header = TRUE, sep = ",")
metaData

#Create the DEseq2 object; assign horse ID in columns w/ metadata in corresponding rows
dds <- DESeqDataSetFromMatrix(
  countData = genecounts, 
  colData = metaData,
  design = ~ timepoint + horse)

#view object
dds

#Filter out genes w/ counts < 10 reads and that are expressed in less than 9 horses, ideally ~50% of sample
dds <- estimateSizeFactors(dds)
keep <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 9

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

#reorder results based on padjust
res <- res[order(res$padj),]
head(res)

#Extract the results and save as a .csv
sig_res <- data.frame(res) %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.05) %>% 
  arrange(padj)

write.csv(
  sig_res, quote = FALSE, row.names = FALSE,
  file = paste0( ResultFile, "SDEGs.csv")
)

#we can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
par(mfrow=c(2,3))

plotCounts(dds, gene="MHCB3", intgroup="timepoint")
plotCounts(dds, gene="ENSECAG00000003925", intgroup="timepoint")
plotCounts(dds, gene="ASS1 ", intgroup="timepoint")
plotCounts(dds, gene="ENSECAG00000033152", intgroup="timepoint")
plotCounts(dds, gene="P3H3", intgroup="timepoint")
plotCounts(dds, gene="TUBGCP5", intgroup="timepoint")

#Volcano Plot Generate
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#PCA Plot #First we need to transform the raw count data
#vst function will perform variance stabilizing transformation

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="timepoint") 

#Sig DEG heat map
mat <- counts(dds, normalized = TRUE)[sig_res$gene, ]
mat_scaled <- t(scale(t(mat)))
M <- (1 - cor(t(mat_scaled), method = "pearson")) / 2
hc <- hclust(as.dist(M), method = "complete")

#Generate heatmap
ht <- Heatmap(
  name = "Scaled\nexpression",
  mat_scaled, 
  col = rev(rainbow(10)),
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
  separate(rowname, remove = F, sep = "_", into = c("timepoint", "horse")) %>%
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
