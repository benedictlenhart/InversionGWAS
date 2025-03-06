### libraries
library(GMMAT)
library(data.table)
library(tidyverse)
library(foreach)
#library(ggplot2)



args = commandArgs(trailingOnly=TRUE)
#create job id from array
id = as.numeric(args[1])
#id = 100
#id = 1#use this when testing to see if script functions correctly


# saveRDS(task.id, paste0(INPUT, "newphenotaskidchart"))
INPUT = "/scratch/bal7cg/fall2023objects/"#path of working directory
ingds = paste0(INPUT,"dgrp.gds")
#load grm
grm = readRDS(paste0(INPUT, "GRM.Aug"))
#load in our phenotype file 
phenotable = readRDS(paste0(INPUT, "sepphenolongform"))
#phenotable = readRDS(paste0(INPUT, "sepphenowideform"))
#filter down to phenotype of interet
phenolist = unique(phenotable$fullpheno)

#phenolist = colnames(phenotable)[-1]
#find gender
# gender = tstrsplit(phenolist[id], "_")
# gender = unlist(tail(gender, n = 1))
#check which wol status to use

dt = phenotable[fullpheno == phenolist[id]]
dt = na.omit(dt)#remove na's 
#fix ral ids
dt = dt %>% 
  dplyr::rename("ral_id" = DGRP) %>% 
  mutate(ral_id = gsub("DGRP", "line", ral_id))
check = dt


#load in wolbachia status data, merge with phenotype data
wol = fread(paste0(INPUT, "female.nStrain20plus.10pca.wolbachia.txt"))
wol = wol[,c("FID","wolbachia"),with=F]

setnames(wol,"FID","ral_id")

setkey(dt,ral_id)
setkey(wol,ral_id)

dt = merge(dt,wol)


#re organize
colnames(dt)[3] <- "avg"
#

#GMMAT creates a null model that uses wolbachia status as fixed effect on phenotype, and the identity matrix as a (nonexistent) random effect. family refers to the method used to calculate the model
#fit a GLMM to the data
modelqtl <- glmmkin(fixed =avg ~ wolbachia, data = dt, kins = grm, id = "ral_id",family = gaussian(link = "identity"))



#run GMMAT
#create  path for output gmmat file to be saved to- using the phenotype variable allows it to be saved to the corresponding directory
outputs = paste(INPUT , 
                "OctFullGRMGWAS/", 
                phenolist[id], ".txt",
                sep = "")


#the glmm.score function scores the impact of each snp agains the null model. minor allele frequency and missing data cutoff can be specified to reduce unwanted computation. compututation can also be sped up by specifying the number of snps calculated per batch, and how many cores to use
glmm.score(modelqtl, infile = ingds, outfile = outputs, MAF.range = c(0.05,0.95), miss.cutoff = 0.15, 
           nperbatch = 400, ncores = 5)


print ("done")
#replace file with rds file to save space
x = fread(paste0("/scratch/bal7cg/fall2023objects/OctFullGRMGWAS/",  phenolist[id], ".txt"))
saveRDS(x, paste(INPUT , 
                 "OctFullGRMGWAS/", 
                 phenolist[id],
                 sep = ""))
file.remove(outputs)#remove txt file
print("finished")

