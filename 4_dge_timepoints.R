#!/usr/bin/Rscript

## Load libraries
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
outName <- "dge"

### Run 
## (1) Run DEseq2 using LRT
## (2) Complete hclus to identify DEG trends
## (3) Attempt correlation analysis


### (1) Run DEseq2 using LRT
#Load in the matrix
mixture <- read.csv("../input_bulk/counts.csv", row.names = 1)
#Make the sample metadata
meta <- data.frame(sample = colnames(mixture)) %>%
  separate(sample, remove = F, sep = "_", into = c("horse", "timepoint"))
#Create the DEseq2 object
dds <- DESeqDataSetFromMatrix(
  mixture, 
  colData = meta,
  design = ~ timepoint + horse
)
#Transform and plot the data with PCA
rld <- varianceStabilizingTransformation(dds, blind=TRUE)
p <- plotPCA(rld, intgroup = "timepoint")
ggsave(paste0("../output/", outName, "/pca.png"), width = 7, height = 7)
p <- plotPCA(rld, intgroup = "horse")
ggsave(paste0("../output/", outName, "/pca2.png"), width = 7, height = 7)
#QC heatmap
rld <- rlog(dds, blind = FALSE)
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
p <- pheatmap::pheatmap(rld_cor)
ggsave(p, file = paste0("../output/", outName, "/pheatmap.png"))
#Complete DE using LRT
dds <- DESeq(dds, test = "LRT", reduced = ~horse)
res <- results(dds)
#Extract the results and save as a .csv
sig_res <- data.frame(res) %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.05) %>% 
  arrange(padj)
write.csv(
  sig_res, quote = FALSE, row.names = FALSE,
  file = paste0("../output/", outName, "/all_genes.csv")
)

### (2) Complete hclus to identify DEG trends
mat <- counts(dds, normalized = TRUE)[sig_res$gene, ]
mat_scaled <- t(scale(t(mat)))
M <- (1 - cor(t(mat_scaled), method = "pearson")) / 2
hc <- hclust(as.dist(M), method = "complete")
ht <- Heatmap(
  name = "Scaled\nexpression",
  mat_scaled, 
  cluster_columns = F, 
  show_row_names = F,
  cluster_rows = F,
  row_split = cutree(hc, k = 6)
)
png(file = paste0("../output/", outName, "/DE_heat.png"), width = 4000, height = 4000, res = 400)
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
#Create the plot
ggplot(data = avg_exp, aes(x = timepoint, y = MEAN)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(clade ~ ., nrow = 2)
ggsave(file = paste0("../output/", outName, "/gex_over_time.png"))
#Print the genes from the clade that decrease with time
print(unique(avg_exp[avg_exp$clade == 2, ]$gene))
#Save the clustering resutlts to complete enrichment scoring w/ scrna data
write.csv(
  avg_exp, quote = FALSE, row.names = FALSE,
  file = paste0("../output/", outName, "/deg_clades.csv")
)

### (3) Compare to synovial fluid -- minimal overlap
res_sf <- read.csv("../../sf/output/allCells_clean/lrt/All cells/allCells_clean_cluster_All cells_all_genes.csv")
mat <- mat[intersect(res_sf$gene, sig_res$gene), ]


### (4) Attempt correlation analysis
lapply(rownames(mat), FUN = function(x){
  browser()
  dat <- data.frame(value = mat[x, ]) %>%
    rownames_to_column() %>%
    separate(rowname, remove = F, into = c("horse", "timepoint")) %>%
    mutate(
      timepoint = as.numeric(gsub("Day", "", timepoint))
    )
  mod <- lm(value ~ timepoint, data = dat)
  ggplot(dat, aes(x=timepoint, y=value)) + 
                geom_smooth(se  = F) + 
    geom_point()
  ggsave(file = paste0("../output/", outName, "/lm.png"))

  lapply(unique(dat$horse), function(y){
    browser()
    dat %>%
      filter(horse == y) %>%
      pivot_wider(names_from = timepoint, values_from = value)
  })
    
})
