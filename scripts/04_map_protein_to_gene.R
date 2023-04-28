## ---------------------------
##
## Script name: 04_map_protein_to_gene
##
## Purpose of script: Maps the protein ids to gene ids
##
## Author: Dr. Veronika Schäpertöns
##
## Date Created: 28.04.2023
##
## Copyright (c) Veronika Schäpertöns, 2023 & Nikolaus Fortelny, 2022
## Email: veronika.schaepertoens@plus.ac.at
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(data.table)


# load data ---------------------------------------------------------------

res <- readRDS(file = "analysis/03_analyse_dge/results_dge.rds")

pmap <- readRDS(file = "analysis/pmap.rds")


# map protein to gene  ----------------------------------------------------

#copy the original data to keep it unchanged and for checking the dimensions after mapping
resOrig <- copy(res)
res <- resOrig

# all are sp entries
stopifnot(all(grepl("^sp\\|", res$rn)))

# Each UniProt ID matches only one gene
stopifnot(all(pmap[,length(unique(Gene)), by = "Entry"][order(V1)]$V1) == 1)
stopifnot(all(pmap[,length(unique(Organism)), by = "Entry"][order(V1)]$V1) == 1)

# save Uniprot ID in a new column
res[,Uniprot := gsub(pattern = "sp\\|(.+?)\\|.+$", replacement = "\\1", x = rn)]

# match res and pmap by Uniprot ID
res <- merge(res, 
             unique(pmap[,c('Entry', "Gene", "Organism"),with = F]), 
             all.x = TRUE, 
             by.x = "Uniprot", 
             by.y = "Entry")

# check that matched res has the same entries as the original
stopifnot(all(res$res %in% resOrig$res))
stopifnot(nrow(res) == nrow(resOrig))
res$res <- NULL

saveRDS(res, "analysis/03_analyse_dge/results_dge_gene.rds")
















