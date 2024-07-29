## ---------------------------
##
## Script name: 02_explore_data
##
## Purpose of script: Explores data with PCA and Pearson's correlation, corrects for batch-effect
##
## Author: Dr. Veronika Schäpertöns & Prof. Nikolaus Fortelny
##
## Date Created: 26.04.2023
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

library(pheatmap)
library(limma)
library(ggplot2)

# load data ---------------------------------------------------------------

load("analysis/all_data.RData")

# plot density and boxplot ------------------------------------------------

plot(density(data.matrix))

boxplot(data.matrix,
        las = 2,
        ylab = "log2 intensity")

# plot PCA ----------------------------------------------------------------

mds <- plotMDS(x = data.matrix,
              col = c(rep("green",9), rep("red",9), rep("blue", 9)),
              labels = data_meta$experiment, 
              gene.selection = "common",
              var.explained = TRUE)

var_explained <- as.data.frame(mds$var.explained[1:7]*100) #first 7 components
colnames(var_explained) <- c("variance")
ggplot(var_explained, aes(x = rownames(var_explained), y = variance)) +
  geom_col() +
  geom_text(aes(label = round(variance, digits = 2)),vjust = -0.2) +
  xlab("principal component number") +
  ylab("variance (%)")

# plot heatmap of Pearson correlations ------------------------------------

pheatmap(cor(data.matrix,method = "spearman"))


# batch-effect correction -------------------------------------------------
data.matrix.batch <- removeBatchEffect(x = data.matrix,
                                       batch = data_meta$experiment)

save(data.matrix, data_meta, data.matrix.batch, file = "analysis/all_data.RData")

# batch-corrected plot PCA ------------------------------------------------

mds <- plotMDS(x = data.matrix.batch,
               col = c(rep("green",9), rep("red",9), rep("blue", 9)),
               labels = data_meta$sample, 
               gene.selection = "common",
               var.explained = TRUE)

var_explained <- as.data.frame(mds$var.explained[1:7]*100) #first 7 components
colnames(var_explained) <- c("variance")
ggplot(var_explained, aes(x = rownames(var_explained), y = variance)) +
  geom_col() +
  geom_text(aes(label = round(variance, digits = 2)),vjust = -0.2) +
  xlab("principal component number") +
  ylab("variance (%)")

# batch-corrected plot heatmap of Pearson correlations --------------------

pheatmap(cor(data.matrix.batch,method = "spearman"))

