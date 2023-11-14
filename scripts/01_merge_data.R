## ---------------------------
##
## Script name: 01_merge_data
##
## Purpose of script: Loads and merges all datasets, loads meta data, saves merged dataset and meta data 
##
## Author: Dr. Veronika Schäpertöns
##
## Date Created: 2023-03-31
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


# load libraries ----------------------------------------------------------

library(tidyverse)


# import datasets ---------------------------------------------------------

data_exp1 <- read_csv("data/Experiment1_Protein_Groups_log2_Normalized_substract_Median_only_valid_values.csv",
                      col_select = c(1:9,30))

data_exp2 <- read_csv("data/Experiment2_Protein_Groups_log2_Normalized_substract_Median_only_valid_values.csv",
                      col_select = c(1:9,30))

data_exp3 <- read_csv("data/Experiment3_new_Protein_Groups_log2_Normalized_substract_Median_only_valid_values.csv",
                      col_select = c(1:9,30))

data_meta <- read_csv("data/pData.csv")

# merge datasets by protein ids -------------------------------------------

data.matrix <- data_exp3 %>% 
  left_join(data_exp1, by = "protein_id") %>% 
  left_join(data_exp2, by = "protein_id") %>%
  na.omit() %>%
  tibble::column_to_rownames("protein_id") %>%
  as.matrix()

dim(data.matrix)

colnames(data.matrix) <- data_meta$sample #rename columns of data.matrix

# save data ---------------------------------------------------------------

save(data.matrix, data_meta, file = "analysis/all_data.RData")




















