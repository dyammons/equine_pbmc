#!/usr/bin/Rscript

## Load libraries
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")
# library(Seurat)
# library(tidyverse)
# library(DESeq2)
# library(ComplexHeatmap)
outName <- "dge"

### This script is split into x sections. Below are detail of each block.
## (1) Run DEseq2 using LRT
## (2) Complete hclus to identify DEG trends
## (3) Run GSEA of the LRT clades
## (4) An aside: compare to ZW analysis
## (5) Reverse CIBERSORT analysis using LRT gene signatures
## (6) Attempt to compare to synovial fluid -- minimal overlap
## (7) Attempt correlation analysis


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
#Complete hclus
mat <- counts(dds, normalized = TRUE)[sig_res$gene, ]
mat_scaled <- t(scale(t(mat)))
M <- (1 - cor(t(mat_scaled), method = "pearson")) / 2
hc <- hclust(as.dist(M), method = "complete")
#Generate heatmap
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
  )
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


### (3) Run GSEA of the LRT clades
avg_exp <- avg_exp[! duplicated(avg_exp$gene), ]
degs <- split(avg_exp$gene, avg_exp$clade)
#Use clusterProfiler
gene_sets <- as.data.frame(
  msigdbr(species = "horse", category = "C5", subcategory = NULL)
) %>% 
  distinct(gs_name, gene_symbol)
res <- lapply(degs, enricher, TERM2GENE = gene_sets)
res <- lapply(res, as.data.frame)
res <- bind_rows(res, .id = 'clade')
res$ID <- gsub("REACTOME_", "", res$ID)

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
png(file = paste0("../output/", outName, "/DE_heat.png"), width = 4500, height = 4000, res = 400)
par(mfcol = c(1, 1))
draw(ht)
dev.off()


### (4) An aside: compare to ZW analysis
#Load in the DE analysis for d14 vs d0 that I had access to
deg.df <- read.csv("../output/zw_analysis/d14_v_d0.csv", row.names = 1)
deg.df <- deg.df[deg.df$padj < 0.05 & abs(deg.df$log2FoldChange) > 2, ]
deg.df$direction <- ifelse(deg.df$log2FoldChange > 2, "UP", "DOWN")
degs <- split(rownames(deg.df), deg.df$direction)
res <- lapply(degs, enricher, TERM2GENE = gene_sets)
res <- lapply(res, as.data.frame)


### (5) Attempt reverse CIBERSORT analysis using LRT gene signatures
## Load the PBMC scrna reference dataset as done in `1_prep_cibersort.R
seu.obj <- readRDS(
    "../../external_data/2023_09_11_pbmc_data_res0.4_dims50_dist0.5_neigh50_S3.rds"
)
#Remove unassigned cells and poorly defined cell types
seu.obj <- subset(
    seu.obj, invert = T, cells = colnames(seu.obj)[is.na(seu.obj$celltype.l2)]
)
seu.obj <- subset(
    seu.obj, invert = T, 
    subset = celltype.l2 %in% c(
        "T CD4- CD8- CD200+",  "B excluded",
        "B proliferating", "T proliferating"
    )
)
#Reduce the number of cell types to improve deconvolution performance
newNames <- structure(
    c(
        "B cell", "CD4 T cell", "CD4 T cell", "CD8 T cell", "CD8 T cell", 
        "CD4 T cell", "CD8 T cell", "B cell", "Monocyte", "B cell",
        "CD8 T cell", "Dendritic cell", "T gd", "Basophil", "Dendritic cell",
        "Dendritic cell", "B cell", "Monocyte", "Neutrophil"
     ), 
    names = 
    c(
        "B T-bet+", "T CD4+ non-naïve", "T CD4+ naïve", "PRF1+ non-annotated",
        "T CD8+ memory", "T CD4+ cytotoxic", "T CD8+ naïve", "B naïve", 
        "classical mono", "B memory", "NK", "cDC2", "T gd", "Basophil", "cDC1",  
        "plasmacytoid DC", "B antibody secreting", "non-classical mono", 
        "Neutrophil"
    )
)
seu.obj <- RenameIdents(seu.obj, newNames)
seu.obj$celltype_cibersort <- Idents(seu.obj)
#Load data and remove duplicate genes
df_clades <- read.csv("../output/dge/deg_clades.csv")
df_clades <- df_clades[! duplicated(df_clades$gene), ]
modulez <- split(df_clades$gene, df_clades$clade)
#Complete module scoring
seu.obj <- AddModuleScore(seu.obj, features = modulez, name = "_score")
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)
#Plot the results
features <- names(modulez)
ecScores <- majorDot(
  seu.obj = seu.obj, groupBy = "celltype_cibersort", 
  scale = T, features = features
) + 
  theme(
    axis.title.y = element_blank(),
    legend.direction = "vertical",
    legend.position = "right"
  ) + 
  labs(x = "Clade") + 
  guides(
    color = guide_colorbar(title = 'Scaled\nenrichment\nscore'),
    size = guide_legend(nrow = 3, byrow = F, title = 'Percent\nenriched')
  )
ggsave(paste0("../output/", outName, "/", outName, "_dots_celltypes.png"), width = 5, height = 4)


### (6) Compare to synovial fluid -- minimal overlap
res_sf <- read.csv("../../sf/output/allCells_clean/lrt/All cells/allCells_clean_cluster_All cells_all_genes.csv")
mat <- mat[intersect(res_sf$gene, sig_res$gene), ]


### (7) Attempt correlation analysis
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
