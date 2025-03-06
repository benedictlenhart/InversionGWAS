library(RhpcBLASctl)
blas_set_num_threads(1)
### libraries
library(GMMAT)
library(data.table)
library(tidyverse)
library(foreach)
#library(ggplot2)



args = commandArgs(trailingOnly=TRUE)
#create job id from array
id = as.numeric(args[1])
#id = 11
#id = 1#use this when testing to see if script functions correctly


# saveRDS(task.id, paste0(INPUT, "newphenotaskidchart"))
INPUT = "/project/berglandlab/Adam/fall2023objects/"#path of working directory
ingds = paste0(INPUT,"dgrp.gds")
#load grm
#grm = readRDS(paste0(INPUT, "LDGRM"))
#load in our phenotype file 

phenotable = readRDS(paste0(INPUT, "pca.long.imputed4"))
#remake pca data to fit older format
phenotable = phenotable %>% 
  mutate(fullpheno = paste(Inversions, Dim, sep = "-")) %>% 
  rename("avg" = "loadings") %>% 
  as.data.table(.)
#filter down to phenotype of interet
phenolist = unique(phenotable$fullpheno)

#check which wol status to use

dt = phenotable[fullpheno == phenolist[id]]
dt = na.omit(dt)#remove na's 
#fix ral ids
dt = dt %>% 
  dplyr::rename("ral_id" = DGRP) %>% 
  mutate(ral_id = gsub("DGRP", "line", ral_id))
check = dt

# 
# #load in wolbachia status data, merge with phenotype data
# wol = fread(paste0(INPUT, "female.nStrain20plus.10pca.wolbachia.txt"))
# wol = wol[,c("FID","wolbachia"),with=F]
# 
# setnames(wol,"FID","ral_id")
# 
# setkey(dt,ral_id)
# setkey(wol,ral_id)
# 
# dt = merge(dt,wol)

#make sure all dgrp lines are in the grm

#

chroms= c("2L", "2R", "3L", "3R", "X")
chrom.out = foreach( f = chroms) %do% {
  # f = "2L"
  #load in the grm for that chromosome
  grm = readRDS(paste0(INPUT, "No", f, "GRM"))
  #assign gds filename based on chromosome 
  ingds = paste0(INPUT,f,".gds")
  dt = dt %>% 
    filter(ral_id %in% colnames(grm))
  succeed = dt$ral_id %in% colnames(grm)
  #fit a GLMM to the data
  
  modelqtl <- glmmkin(fixed =avg ~ 1, data = dt[!is.na(avg)], kins = grm, id = "ral_id",family = gaussian(link = "identity"))
  #create  path for output gmmat file to be saved to- using the phenotype variable allows it to be saved to the corresponding directory
  
  outputs = paste(INPUT , 
                  "loco.pca.gwas4/", 
                  phenolist[id], "#",f,"chrom.txt",
                  sep = "")
  # outputs = paste(INPUT , 
  #                 "loco.pca.gwas4/", 
  #                 "test.txt",
  #                 sep = "")
    
  
  glmm.score(modelqtl, infile = ingds, outfile = outputs, MAF.range = c(0.05,0.95), miss.cutoff = 0.15, 
             nperbatch = 400, ncores = 5)
  x = fread(outputs)
  saveRDS(x, paste(INPUT , 
                   "loco.pca.gwas4/", 
                   phenolist[id],"#",f,
                   sep = ""))
  file.remove(outputs)#remove txt file
  print(f)
  print("finished")
  
}



