.## ---------------------------
##
## Script name: 07_plot_gsea_results
##
## Purpose of script: Plots NES and padj of GSEA results
##
## Author: Dr. Veronika Schäpertöns & Prof. Nikolaus Fortelny
##
## Date Created: 01.09.2023
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
library(tidyverse)
library(data.table)
library(RColorBrewer)


# load gsea results -------------------------------------------------------


gsea.res <- read_tsv("analysis/06_analyse_gsea/GSEA.tsv") %>%
  data.table()



dim(gsea.res[padj < 0.05])

gsea.res[padj < 0.05]

data_to_plot <- gsea.res[padj < 0.05][,-"leadingEdge"][grp == "drfae_vs_hpyl"][order(padj)]

write.table(data_to_plot, "analysis/06_analyse_gsea/gsea_drfae_hpyl.tsv")

data_to_plot <- read_delim("analysis/06_analyse_gsea/gsea_drfae_hpyl_short.tsv",
                           delim = " ")


# bar graph ---------------------------------------------------------------

ggplot(data = data_to_plot) +
  geom_col(
    mapping = aes(x = NES,
                  y = pathway,
                  fill = -log10(padj))
  ) +
  scale_fill_distiller(type = "seq", palette = "YlGn", direction = 1) 


# dot plot ----------------------------------------------------------------

#svg(filename = "figures/gsea_drfae_wt.svg")

ggplot(data_to_plot, aes(x = "drfae_vs_wt", y = pathway, size = -log10(padj))) +
  geom_point(aes(color = NES)) +
  scale_color_distiller(type = "div", palette = "PRGn") +
  xlab("comparison (drfae vs wt)") +
  ylab("gene sets") 

#dev.off()

