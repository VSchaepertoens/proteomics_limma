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
library(paletteer)


# load gsea results -------------------------------------------------------


gsea.res <- read_tsv("analysis/GSEA.tsv") %>%
  data.table()



dim(gsea.res[padj < 0.05])

gsea.res[padj < 0.05]

data_to_plot <- gsea.res[padj < 0.05][,-"leadingEdge"][grp == "drfae_vs_hpyl"][order(padj)]

write.table(data_to_plot, "analysis/06_analyse_gsea/gsea_drfae_hpyl.tsv")

data_to_plot <- read_delim("analysis/06_analyse_gsea/gsea_drfae_hpyl_short.tsv",
                           delim = " ")
 

# bar graph ---------------------------------------------------------------

svg(filename = "figures/bargraph_gsea_drfae_wt_darkRed.svg", width = 2.8, height = 2.6)

ggplot(data = data_to_plot) +
  geom_col(
    mapping = aes(x = NES,
                  y = reorder(pathway, NES),
                  fill = -log10(padj))
  ) +
  scale_fill_paletteer_c("ggthemes::Classic Red-Black", direction = -1) +
  xlab("normalized enrichment score") +
  ylab("gene sets") +
  scale_fill_gradient(low = "black", high = rgb(0.62,0,0), name = expression(-log[10](p[adj]))) +
  theme_bw() +
  theme(text = element_text(size = 6, family = "sans"),
        legend.key.height = unit(4, 'mm'),
        legend.key.width = unit(2, 'mm'),
        legend.title = element_text(angle = 90),
        legend.box.spacing = unit(1, "pt")) # The spacing between the plotting area and the legend box (unit))
  
#RGB (160, 0, 0) --> rgb(0.62, 0, 0)

dev.off()  
# dot plot ----------------------------------------------------------------

#svg(filename = "figures/gsea_drfae_wt.svg")
svg(filename = "figures/dotplot_gsea_drfae_wt.svg", width = 2.8, height = 2.6)


ggplot(data_to_plot, aes(x = "drfae_vs_wt", y = reorder(pathway, NES), size = -log10(padj))) +
  geom_point(aes(color = NES)) +
  scale_color_distiller(type = "div", palette = "PRGn") +
  xlab("comparison (drfae vs wt)") +
  ylab("gene sets") +
  scale_size(range = c(0, 3)) +
  theme_bw() +
  theme(text = element_text(size = 6, family = "sans"),
        legend.key.height = unit(4, 'mm'),
        legend.key.width = unit(2, 'mm'),
        #legend.title = element_text(angle = 90),
        legend.box.spacing = unit(1, "pt")) # The spacing between the plotting area and the legend box (unit))


#dev.off()
dev.off()  

