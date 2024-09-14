#!/usr/bin/Rscript

## Load libraries
library(Seurat)


### This script is split into 2 sections. Below are detail of each block.
###
### Purpose: extract the cell type expression from an equine PBMC dataset
### (obtained from Patel et al.) and prepare out bulk seq count matrix to be
### loaded into CIBERSORTx to estimate how cell fractions change
###
### Sections:
## (1) Load the scrna reference, update cell type levels, and export matrix
## (2) Load the bulk seq count matrix and format to use in CIBERSORT


### (1) Load the scrna reference, update cell type levels, and export matrix
## Load the PBMC scrna reference dataset and alter annotation levels
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
#Set the variable to use for cibersort to avoid hardcoding later
ct_level <- "celltype_cibersort"
## Export the data as a matrix for import into CIBERSORTx
#Down sample to make data easier to work with (keep equal number of cells)
Idents(seu.obj) <- ct_level
seu.obj <- subset(
    x = seu.obj, downsample = min(table(seu.obj@meta.data[[ct_level]]))
)
message(paste0(
    "Each cell type was downsampled to ",
    min(table(seu.obj@meta.data[[ct_level]])),
    " cells.\nConsider this value when selecting replicates ",
    "to build the deconvolution reference."
))
#Get the data - this can use a lot of RAM if there are still a lot of cells
#Here we get raw counts - CIBERSORTx should convert this to CPM
cnt_mat <- FetchData(
    seu.obj, vars = rownames(seu.obj), assay = "RNA", layer = "counts"
)
#Transpose the data to get it in genes (rows) X cells (columns)
cnt_mat <- t(cnt_mat)
#Replace the cell barcodes with the cell type name
colnames(cnt_mat) <- as.character(seu.obj@meta.data[[ct_level]])
cnt_mat <- cbind(rownames(cnt_mat), cnt_mat)
colnames(cnt_mat)[1] <- "GeneSymbol"
cnt_mat <- as.data.frame(cnt_mat)
#Save the matrix
write.table(
    cnt_mat, quote = FALSE, row.names = FALSE, sep = "\t",
    file = paste0("../output/cibersort_data/scrna_", gsub("\\.", "_", ct_level), ".txt")
)

### (2) Load the bulk seq count matrix and format to use in CIBERSORT
#Read input counts data, then write as tsv
mixture <- read.csv("../input_bulk/counts.csv")
write.table(
    mixture, quote = FALSE, row.names = FALSE, sep = "\t",
    file = "../output/cibersort_data/mixture.txt"
)
