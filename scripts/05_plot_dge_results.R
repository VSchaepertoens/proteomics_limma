## ---------------------------
##
## Script name: 05_plot_dge_results
##
## Purpose of script: Plots heatmaps, pval histogram, volcano plots
##
## Author: Dr. Veronika Schäpertöns & Prof. Nikolaus Fortelny
##
## Date Created: 28.04.2023
##
## Copyright (c) Veronika Schäpertöns, 2023 
## Email: veronika.schaepertoens@plus.ac.at
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(ggplot2)
library(ggrepel)
library(tidyverse)
library(pheatmap)

# load data ---------------------------------------------------------------

load(file = "analysis/results_dge.RData")

load(file = "analysis/all_data.RData")

# plot results ------------------------------------------------------------

# vulcano plots
ggplot(res, aes(x = logFC, y = -log10(P.Value))) +
  geom_point() + 
  facet_wrap(~coef)

ggsave("figures/figure_4b.pdf")


# p-value histograms
ggplot(res, aes(x = P.Value)) +
  geom_histogram() + 
  facet_wrap(~coef)

## compare drfae vs hpyl
# extract logFC values of drfae vs.untreated with hpyl vs. untreated
pDT <- merge(suffixes = c("_drfae", "_hpyl"), by = "rn",
             res[coef == "drfae"][, c("rn", "logFC"), with = F],
             res[coef == "hpyl"][, c("rn", "logFC"), with = F]
)

# extract adj.P.Val values of drfae vs. hpyl 
pDT <- merge(pDT, res[coef == "drfae_vs_hpyl"][,c("rn", "adj.P.Val")], by = "rn")

ggplot(pDT[adj.P.Val < 0.05], aes(x = logFC_hpyl, y = logFC_drfae)) + 
  #geom_hex() + 
  geom_text_repel(data = pDT[adj.P.Val < 0.05][!grepl("\\;", rn)], aes(label = gsub("^.+\\|(.+?)_HUMAN$", "\\1", rn))) + 
  geom_point(data = pDT[adj.P.Val < 0.05], shape = 1) + 
  scale_fill_gradient(low = "grey", high = "#1f78b4") + 
  theme_bw(12) + 
  geom_abline()


# plot heatmaps of log2 values --------------------------------------------
## function to plot top genes
plot_top_genes <- function(mat, 
                           comparison,
                           top_genes = FALSE,
                           no_genes = 20,
                           show_rownames = FALSE){
  
  ann_df <- data.frame(logFC = res_matched[coef == comparison][adj.P.Val < 0.05]$logFC,
                       rn = res_matched[coef == comparison][adj.P.Val < 0.05]$rn,
                       gene = res_matched[coef == comparison][adj.P.Val < 0.05]$Gene)
  
  if (top_genes) {
    ann_df <- ann_df[order(abs(ann_df$logFC), decreasing = TRUE),]  
    ann.row <- with(ann_df[1:no_genes,], 
                  data.frame(row.names = gene,log2FC = ifelse(logFC > 0, "up","down")))
  
    pheatmap(mat[ann_df$gene[1:no_genes], ], 
             annotation_colors = list(log2FC = c(down = "blue",up = "red")),
             annotation_row = ann.row,
             cluster_cols = T,
             show_rownames = show_rownames,
             scale = "row")
  } else {
    ann.row <- with(ann_df,
                    data.frame(row.names = gene, logFC = ifelse(logFC > 0, "up", "down")))
    pheatmap(mat[ann_df$gene,],
             annotation_colors = list(logFC = c(down = "blue",up = "red")),
             annotation_row = ann.row,
             cluster_cols = TRUE,
             show_rownames = FALSE,
             scale = "row") 
    }
}

# heatmap of drfae vs untreated, all significant genes, all samples
plot_top_genes(data.matrix.batch.matched, "drfae")

# heatmap of hpyl vs untreated, all significant genes, all samples
plot_top_genes(data.matrix.batch.matched, "hpyl")

# heatmap of drfae vs hpyl, all significant genes, all samples
plot_top_genes(data.matrix.batch.matched, "drfae_vs_hpyl")


# subset only treated batch-corrected samples -----------------------------

data.matrix.batch.matched.treated <- data.matrix.batch.matched[, -grep("un", 
                                                                       colnames(data.matrix.batch.matched))]

# heatmap of drfae vs hpyl, top 20 significant genes, treated samples
plot_top_genes(mat = data.matrix.batch.matched.treated, 
               comparison = "drfae_vs_hpyl",
               top_genes = TRUE, 
               no_genes = 20, 
               show_rownames = TRUE)

# heatmap of drfae vs hpyl, top 40 significant genes, treated samples
plot_top_genes(data.matrix.batch.matched.treated, 
               "drfae_vs_hpyl", 
               TRUE,
               40, 
               TRUE)

# heatmap of drfae vs hpyl, top 60 significant genes, treated samples
plot_top_genes(data.matrix.batch.matched.treated,
               "drfae_vs_hpyl", 
               TRUE,
               60, 
               TRUE)

# heatmap of drfae vs hpyl, all significant genes, treated samples
pdf("figures/figure_4a.pdf")
plot_top_genes(data.matrix.batch.matched.treated, 
               "drfae_vs_hpyl", 
               TRUE,
               427, 
               FALSE) 

dev.off()


