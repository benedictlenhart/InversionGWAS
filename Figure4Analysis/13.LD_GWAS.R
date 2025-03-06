### libraries
library(GMMAT)
library(data.table)
library(tidyverse)
library(foreach)


args = commandArgs(trailingOnly=TRUE)
#create job id from array
id = args[1]
job = 1

#id = 1#use this when testing to see if script functions correctly


# saveRDS(task.id, paste0(INPUT, "newphenotaskidchart"))
INPUT = "/scratch/bal7cg/summer2023objects/"#path of working directory



#load in the phenotype table

perms.tables = readRDS(paste0(INPUT,"wideform.fixed.phenotable.RDS"))

#get colnames of phenotypes
phenonames = colnames(perms.tables)[2:length(colnames(perms.tables))]
phenotype = phenonames[job + 1]#identify the phenotype assaye dhere
#define columns we want to keep
cols = c("ral_id",phenotype)
dt = perms.tables[,cols,with=F]



#load in wolbachia status data, merge with phenotype data
wol = fread(paste0(INPUT,"female.nStrain20plus.10pca.wolbachia.txt"))
wol = wol[,c("FID","wolbachia"),with=F]
#remove line_ from wol
dt = dt %>% 
  mutate(ral_id = paste("line", ral_id, sep = "_"))
setnames(wol,"FID","ral_id")

setkey(dt,ral_id)
setkey(wol,ral_id)

dt = merge(dt,wol)

#re organize
colnames(dt)[2] <- "avg"

grm = readRDS(paste0(INPUT, "LDGRM"))
#assign gds filename based on chromosome 
ingds = paste0("dgrp.gds")

#fit a GLMM to the data

modelqtl <- glmmkin(fixed =avg ~ wolbachia, data = dt[!is.na(avg)], kins = grm, id = "ral_id",family = gaussian(link = "identity"))
#create  path for output gmmat file to be saved to- using the phenotype variable allows it to be saved to the corresponding directory
outputs = paste0(INPUT , "ldprunegwas/", 
                 phenotype,".txt")
glmm.score(modelqtl, infile = ingds, outfile = outputs, MAF.range = c(0.05,0.95), miss.cutoff = 0.15, 
           nperbatch = 400, ncores = 5)

print ("done")


