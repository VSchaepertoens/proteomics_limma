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
library(pheatmap)

# load data ---------------------------------------------------------------

data.matrix.batch <- readRDS("analysis/02_explore_data/batch_corrected_all_data.rds")

res <- readRDS("analysis/03_analyse_dge/results_dge_gene.rds")

load(file = "analysis/02_explore_data/matched_treated_batch_corrected_all_data.RData")


# plot results ------------------------------------------------------------

# vulcano plots
ggplot(res, aes(x = logFC, y = -log10(P.Value))) +
  geom_point() + 
  facet_wrap(~coef)

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

# heatmap of drfae vs untreated, all significant genes, all samples
plot_genes <- function(mat,
                       comparison,
                       detailed_annotation = FALSE,
                       top_n = Inf) {
  selected_genes <- 
    res[coef == "drfae_vs_hpyl"] %>%
    as_tibble() %>% 
    mutate(abs_logFC = abs(logFC)) %>% 
    slice_max(n = top_n, abs_logFC) %>% 
    pull(Gene)
  
  if (detailed_annotation) {
    annx <- dcast.data.table(res, rn ~ coef, value.var = "logFC")
    ann.row <- data.frame(row.names = annx$rn, annx[, -"rn", with = FALSE])
  } else {
    ann.row <- with(
      res[coef == comparison],
      data.frame(row.names = rn, logFC = ifelse(logFC > 0, 1, -1))
    )
  }
  
  pheatmap(mat[res[coef == comparison][adj.P.Val < 0.05]$rn,],
           annotation_colors = list(logFC = c("blue","red")),
           annotation_row = ann.row,
           cluster_cols = TRUE,
           show_rownames = FALSE,
           scale = "row")  
}

plot_genes(data.matrix.batch, "drfae")

# heatmap of hpyl vs untreated, all significant genes, all samples
plot_genes(data.matrix.batch, "hpyl")

# heatmap of drfae vs hpyl, all significant genes, all samples
plot_genes(data.matrix.batch, "drfae_vs_hpyl")

# heatmap of drfae vs hpyl, all significant genes, all samples, ANNOTATE by actual logFC values
plot_genes(data.matrix.batch, "drfae_vs_hpyl", TRUE)

# heatmap of drfae vs hpyl, all significant genes, ONLY TREATED samples



#plot only top genes
ann_df <- data.frame(logFC = res[coef == "drfae_vs_hpyl"][adj.P.Val < 0.05]$logFC,
                     rn = res[coef == "drfae_vs_hpyl"][adj.P.Val < 0.05]$rn,
                     gene = res[coef == "drfae_vs_hpyl"][adj.P.Val < 0.05]$Gene)

ann_df <- ann_df[order(abs(ann_df$logFC), decreasing = TRUE),]  

#ann.row <- with(res[coef == "drfae_vs_hpyl"], data.frame(row.names = rn,logFC = ifelse(logFC > 0, 1,-1)))
ann.row <- with(ann_df[1:20,], data.frame(row.names = gene,logFC = ifelse(logFC > 0, 1,-1)))

#top 20
pheatmap(data.matrix.batch.treated.matched[selected_genes, ], 
         annotation_colors = list(logFC = c("blue","red")),
         annotation_row = ann.row,
         cluster_cols = T,
         show_rownames = TRUE,
         scale = "row")


data.matrix.batch[, ...] %>% 
  plot_genes()

plot_genes(data.matrix.batch[, ...], )




