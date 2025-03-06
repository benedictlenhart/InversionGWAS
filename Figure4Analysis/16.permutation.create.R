library(lme4)
library(data.table)
library(tidyverse)
library(plotrix)
library(gmodels)
library(scales)
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Sep_2023_objects/")

#create 100 permutations of the pheno table with dgrp ids randomly shuffled

#wideform data 
wide = readRDS("sepphenowideform")


# Function to shuffle the first column of a data table
shuffle_first_column <- function(dt) {
  dt[,DGRP := sample(DGRP, replace = FALSE)]
}
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Sep_2023_objects/phenoperms/")
# Shuffle the first column 100 times
for (i in 1:100) {
  set.seed(i)
  #i = 1
  dt <- shuffle_first_column(wide)
  # You can perform any further operations with the shuffled data table here
  saveRDS(dt, paste0("phenopermutation", i))
}

# Print the shuffled data table
print(dt)
#make a chart table that keys which perm and which phenotype to run gwas for
# Create ID variables ranging from 1 to 100 and 200
id_var_1 <- 1:100
id_var_2 <- colnames(wide)
id_var_2 <- id_var_2[-1]
# Create all pairwise combinations using expand.grid
pairs_combinations <- expand.grid(ID_var_1 = id_var_1, ID_var_2 = id_var_2)

# Convert the result to a data.table
ref.table <- as.data.table(pairs_combinations)

48800 / 5 #how many array jobs will we need with 5 gwas per job
ref.table$array = rep(1:9760, each = 5)
saveRDS(ref.table, "/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2024_objects/ref.table")
