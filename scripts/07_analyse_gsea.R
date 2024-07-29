## 
##
## Script name: 07_analyse_gsea
##
## Purpose of script: runs fGSEA analysis 
##
## Author: Dr. Veronika Schäpertöns & Prof. Nikolaus Fortelny
##
## Date Created: 23.03.2023
##
## Copyright (c) Veronika Schäpertöns, 2023
## Email: veronika.schaepertoens@plus.ac.at
##
## 
##
## Notes:
##   
##
## 

library(tidyverse)
library(data.table)
library(fgsea)

# load data ---------------------------------------------------------------

load(file = "analysis/results_dge.RData")

# Number of multi-matches
table(grepl("\\;", res_matched$rn))

nrow(res_matched[!is.na(Gene)])

agg.res.enr <- res_matched[!is.na(Gene) & Organism == "Homo sapiens (Human)", .(score = mean(logFC, na.rm = TRUE)), by = c("Gene", "coef")]
nrow(agg.res.enr)

# load genesets (pathways for fgsea) ------------------------------

genesets <- read_tsv("analysis/enrichr_database.tsv")

# enrichment analysis using fgsea -----------------------------------------
set.seed(119)
gsea.res <- data.table()
de.grp <- res$coef[1]

for (de.grp in unique(res_matched$coef)) {
  print(de.grp)
  gsea.res <- rbind(gsea.res, 
                    data.table(fgsea(
                      pathways = with(genesets, split(Gene, paste(DB, Geneset, sep = "__"))),
                      stats = with(agg.res.enr[coef == de.grp], setNames(score, nm = Gene)),
                      nperm = 1e6), 
                      grp = de.grp))
}

gsea.res_orig <- gsea.res
#gsea.res <- gsea.res_orig

gsea.res$leadingEdge <- sapply(
  gsea.res$leadingEdge, 
  function(vec) paste(vec, collapse = ",")
)


# gsea results ------------------------------------------------------------
dim(gsea.res[padj < 0.05])

gsea.res[padj < 0.05]
gsea.res[padj < 0.05][,-"leadingEdge"][order(padj)]

gsea.res[padj < 0.05][,-"leadingEdge"][grp == "hpyl"][order(padj)]
gsea.res[padj < 0.05][,-"leadingEdge"][grp == "drfae"][order(padj)]
gsea.res[padj < 0.05][,-"leadingEdge"][grp == "drfae_vs_hpyl"][order(padj)]


write_tsv(x = gsea.res, file = "analysis/results_gsea.tsv")

