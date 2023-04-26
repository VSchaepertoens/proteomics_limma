## ---------------------------
##
## Script name: 03_analyse_dge
##
## Purpose of script: Analyses differential protein expression using linear models package limma
##
## Author: Dr. Veronika Schäpertöns & Nikolaus Fortelny
##
## Date Created: 26.04.2023
##
## Copyright (c) Veronika Schäpertöns, 2023 & Nikolaus Fortelny 2022
## Email: veronika.schaepertoens@plus.ac.at
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(data.table)
library(limma)

# load data ---------------------------------------------------------------

load("analysis/01_merge_data/all_data.RData")


# design matrix -----------------------------------------------------------

data_meta_refactor <- copy(data_meta)
unique(data_meta_refactor$treatment)
data_meta_refactor$treatment <- factor(x = data_meta_refactor$treatment, 
                                       levels = c("untreated", "hpyl", "drfae"))
unique(data_meta_refactor$treatment)

design.matrix <- model.matrix(~experiment + treatment, data = data_meta_refactor)
row.names(design.matrix) <- data_meta$sample


# model fit ---------------------------------------------------------------

fit <- lmFit(data.matrix, design.matrix)
fit <- eBayes(fit)


# get results -------------------------------------------------------------

coefs <- grep("treatment", colnames(coef(fit)), value = TRUE)
res <- data.table()
for (coefx in coefs) {
  res <- rbind(res, data.table(
    topTable(fit, coef = coefx, number = nrow(data.matrix)), 
    keep.rownames = TRUE,
    coef = gsub("treatment", "", coefx)
  ))
}
res[coef == "drfae"][adj.P.Val < 0.05]
res[coef == "hpyl"][adj.P.Val < 0.05]


# contrast fit ------------------------------------------------------------

coef.all <- colnames(coef(fit))
coef.used <- c("treatmentdrfae", "treatmenthpyl")
contr.use <- data.frame(
  "drfae_vs_hpyl" = c(c(1, -1), rep(0, length(coef.all) - length(coef.used))), 
  row.names = c(coef.used, setdiff(coef.all, coef.used)))

contrasts = contr.use[coef.all,,drop = F]

fit2 <- contrasts.fit(fit, as.matrix(contrasts))
fit2 <- eBayes(fit2)

coefx <- colnames(contr.use)
res <- rbind(res, data.table(
  topTable(fit2, coef = coefx, number = nrow(data.matrix)), 
  keep.rownames = TRUE,
  coef = gsub("treatment", "", coefx)
))

res[,direction := ifelse(logFC > 0, "up", "down")]

res$res <- 1:nrow(res)
resOrig <- copy(res)

# Number of significant hits
res[adj.P.Val < 0.05][,.N, by = c("coef", "direction")]
# Number of tested
res[,.N, by = c("coef", "direction")]

saveRDS(res, file = "analysis/03_analyse_dge/results_dge.rds")


























