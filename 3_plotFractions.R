#!/usr/bin/Rscript

## Load libraries
library(tidyverse)

### This script is split into 2 sections. Below are detail of each block.
###
### Purpose: identify if cibersort batch correction methods look to yeild 
### conistent results then compare estimated cell type percentages across the
### time series using t.test.
###
### Sections:
## (1) Compare CIBERSORTx batch correction approaches
## (2) Compare cell type estimates between conditions

### (1) Compare CIBERSORTx batch correction approaches
#Load in the results
filez <- c(
    "../output/cibersort_output/noBatch/CIBERSORTx_Results.txt",
    "../output/cibersort_output/bBatch/CIBERSORTx_Adjusted.txt",
    "../output/cibersort_output/sBatch/CIBERSORTx_Adjusted.txt"
)
df.list <- lapply(filez, function(x){
    res <- read.table(file = x, sep = "\t", header = TRUE)
    res$method <- strsplit(x, "/")[[1]][4]
    return(res)
})
res <- do.call(rbind, df.list)
#Prep data for plotting with ggplot
res <- pivot_longer(res, cols = colnames(res)[2:(ncol(res) - 4)])
res$value <- round(res$value, 2)
res$name <- gsub("\\.", " ", res$name)
res <- left_join(res, meta[ , c(1,3)], by = c("Mixture" = "orig.ident"))
res$method <- factor(res$method, levels = c("noBatch", "bBatch", "sBatch"))
### Plot heatmap of estimated percentages
ggplot(data = res, aes(name, Mixture, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
        low = "white", high = "red", limit = c(0,1), space = "Lab", 
        name = "Estimated cell\nfraction"
    ) +
    theme_minimal() + 
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, 
        size = 12, hjust = 1)
    ) +
    facet_grid(.~method) +
    coord_fixed() + 
    geom_text(aes(name, Mixture, label = value), color = "black", size = 2) +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(paste0("../output/cibersort_output/cell_fractions.png"), width = 16, height = 5)
### Conclusion: S batch corrected data looks reasonable

### (2) Compare cell type estimates between conditions
#Subset on sBatch estimates and clean dataframe
res <- res[res$method == "sBatch", ] 
res <- separate(
    res, Mixture, into = c("Horse", "Timepoint"), sep = "_", remove = F
)
res$Timepoint <- factor(res$Timepoint, levels =c("Day0", "Day14", "Day28", "Day56", "Day126"))
stat_res <- res %>% 
    pivot_wider(names_from = name, values_from = value)
## Run t.test between Timepoint for each cell type and annotate data
#Run t.test
stat.dat <- lapply(colnames(stat_res)[8:15], function(x){
    inner.dat <- lapply(levels(res$Timepoint)[2:5], function(y){
        stat_res <- stat_res[stat_res$Timepoint %in% c(y, "Day0"), c(colnames(stat_res)[c(1:3)], x)]
        colnames(stat_res)[4] <- "value"
        if (length(which(stat_res$Timepoint == "Day0")) !=
            length(which(stat_res$Timepoint == y))) {
            stat_res <- filter(stat_res, Horse %in% stat_res[stat_res$Timepoint == y, ]$Horse)
        }
        pval <- t.test(
            stat_res[stat_res$Timepoint == "Day0", ]$value,
            stat_res[stat_res$Timepoint == y, ]$value,
            paired = T, alternative = "two.sided"
        )$p.value
        dat <- as.data.frame(matrix(
            c(x, y, pval),
            dimnames = list(1, c("Cell_type", "TimePoint", "Pvalue")),
            ncol = 3, nrow = 1, byrow = TRUE
        ))
        return(dat)
    })
    do.call(rbind, inner.dat)
})
stat.dat <- do.call(rbind, stat.dat)
stat.dat$Pvalue <- as.numeric(stat.dat$Pvalue)
#Annotate the significant differences
pVal <- 0.1
sig.res <- stat.dat
sig.res$y_pos_dot <- rep(unlist(
    lapply(colnames(stat_res)[8:15], function(ct) {
        max(na.omit(res[res$name == ct, ]$value)) + max(res$value) * 0.1
    })
), each = 4)
sig.res$y_pos_line <- sig.res$y_pos_dot - max(res$value) * 0.04
sig.res$row_num <- rep(1:8, each = 4)
sig.res <- sig.res %>%
  mutate(
      x_start = row_num - 0.4,
      x_end = row_num + 0.4,
      sig_stars = case_when(
          Pvalue < 0.001 ~ "***",
          Pvalue > 0.001 & Pvalue < 0.01 ~ "**",
          Pvalue > 0.01 & Pvalue < 0.05 ~ "*",
          Pvalue < 0.001 ~ "n.s"
      )
  )
sig.res <- sig.res[sig.res$Pvalue < pVal, ]
sig.res <- na.omit(sig.res)
#Create plot with signifcance annotation
res$name <- factor(res$name, levels = colnames(stat_res)[8:15])
ggplot(res, aes(x = name, y = value)) +
    geom_bar(
        aes(fill = Timepoint, colour = Timepoint), 
        alpha = 0.5, stat = "summary", fun = mean, 
        position = position_dodge2(padding = 0.2)
    ) + 
    ggnewscale::new_scale_fill() +
    geom_point(
        aes(fill = Timepoint, group = Timepoint), 
        shape = 21, colour = "grey20", size = 1.5, alpha = 0.8,
        position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1),
        show.legend = FALSE
    ) + 
    { 
      if (sum(sig.res$Pvalue < pVal) > 0) { 
          annotate("text", x = sig.res$Cell_type, y = sig.res$y_pos_dot, label = sig.res$sig_stars)
      }
    } + 
    {
      if (sum(sig.res$Pvalue < pVal) > 0) { 
        geom_segment(
          data = sig.res, aes(x = x_start, xend = x_end, y = y_pos_line, yend = y_pos_line)
        )
      }
    } +
    labs(
      x = "Cell type", 
      y = "CIBERSORTx estimated\ncell fraction (%)",
      title = "Cell fraction over time"
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_blank(),
      legend.key = element_rect(fill = 'transparent', colour = NA),
      legend.background = element_rect(fill = 'transparent', colour = NA)
    ) + 
  coord_cartesian(expand = FALSE, clip = 'off') +
  scale_y_continuous(limits = c(0, max(c(res$value, sig.res$y_pos_dot)) * 1.1))
ggsave("../output/cibersort_output/btwn_grps.png", width = 6, height = 4)
