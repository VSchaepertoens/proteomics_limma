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

load(file = "analysis/results_dge.RData")

load(file = "analysis/all_data.RData")

pmap <- readRDS(file = "analysis/pmap.rds")

## res ##
# map protein to gene  ----------------------------------------------------

#copy the original data to keep it unchanged and for checking the dimensions after mapping
resOrig <- copy(res)
#res <- resOrig

# all are sp entries
stopifnot(all(grepl("^sp\\|", res$rn)))

# Each UniProt ID matches only one gene
stopifnot(all(pmap[,length(unique(Gene)), by = "Entry"][order(V1)]$V1) == 1)
stopifnot(all(pmap[,length(unique(Organism)), by = "Entry"][order(V1)]$V1) == 1)

# save Uniprot ID in a new column
#res[!grepl("\\;", rn),Uniprot := gsub("sp\\|(.+)\\|.+$", "\\1", rn)]
res[,Uniprot := gsub(pattern = "sp\\|(.+?)\\|.+$", replacement = "\\1", x = rn)]

# match res and pmap by Uniprot ID
res <- merge(res, 
             unique(pmap[,c('Entry', "Gene", "Organism"),with = FALSE]), 
             all.x = TRUE, 
             by.x = "Uniprot", 
             by.y = "Entry")

# check that matched res has the same entries as the original
stopifnot(all(res$res %in% resOrig$res))
stopifnot(nrow(res) == nrow(resOrig))
res$res <- NULL

res_matched <- res
res <- resOrig

save(res, res_matched, file = "analysis/results_dge.RData")

## data.matrix.batch ##
# map protein to gene  ----------------------------------------------------

data.matrix.batch_orig <- copy(data.matrix.batch)
data.matrix.batch <- cbind(rn = rownames(data.matrix.batch),
                           data.matrix.batch)
data.matrix.batch <- data.table(data.matrix.batch)
data.matrix.batch[, Uniprot := gsub("sp\\|(.+?)\\|.+$", "\\1", rn)]

data.matrix.batch <- merge(data.matrix.batch,
                           unique(pmap[,c('Entry', "Gene", "Organism"),with = F]),
                           all.x = TRUE,
                           by.x = "Uniprot",
                           by.y = "Entry")

# check for duplicates
idx <- data.matrix.batch[duplicated(data.matrix.batch$Gene) == TRUE, which = TRUE]
data.matrix.batch <- data.matrix.batch[-(idx),]

# filter out only numeric columns
data.matrix.batch.matched <- data.matrix.batch[,-c(1,2,30,31)]
data.matrix.batch.matched <- apply(data.matrix.batch.matched,
                                   MARGIN = 2,
                                   FUN = as.numeric)
rownames(data.matrix.batch.matched) <- data.matrix.batch$Gene

save(data.matrix, data_meta, data.matrix.batch, data.matrix.batch.matched, file = "analysis/all_data.RData")
