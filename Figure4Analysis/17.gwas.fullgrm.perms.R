#***************************
#   The following lines need to keep the glmm.scores function from freezing

library(RhpcBLASctl)
blas_set_num_threads(1)
#***************************
### libraries
library(GMMAT)
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel(5)

#library(ggplot2)

print("libraries loaded")

args = commandArgs(trailingOnly=TRUE)
#create job id from array
id = as.numeric(args[1])
#id = 4
#id = 1#use this when testing to see if script functions correctly
numCores <- as.numeric(args[2]) - 1
#numCores = 4
print(id)
print(numCores)

# saveRDS(task.id, paste0(INPUT, "newphenotaskidchart"))
INPUT = "/scratch/bal7cg/fall2023objects/"#path of working directory
ingds = paste0(INPUT,"dgrp.gds")
#load in reftable
reftable = readRDS(paste0(INPUT,"dgrpreftable"))
#load grm
grm = readRDS(paste0(INPUT, "GRM.Aug"))
#load in array table
arrayref = readRDS(paste0(INPUT, "ref.table"))
#filter to this array's section
arrayref = arrayref[array == id]
#create a foreach loop that loads in the appropriate file and runs gwas per permutation and phenotype
out = foreach(f = as.vector(arrayref[[1]])) %do% {
  #f = 16
  #load in our permuted phenotype file based on array
  phenotable = readRDS(paste0(INPUT,"/phenoperms/phenopermutation",f))
  #melt phenotable
  melt = melt(phenotable, id.vars = ("DGRP"), variable.name = "fullpheno", value.name = "mean")
  #phenotype determined by array table
  pheno.id = as.character(arrayref[ID_var_1 == f]$ID_var_2)
  
  
  dt = melt[fullpheno == pheno.id]
  #merge in the values for reftable
  dt = merge(dt, reftable, by = "DGRP" )
  dt = na.omit(dt)#remove na's 
  #fix ral ids
  dt = dt %>% 
    dplyr::rename("ral_id" = DGRP) %>% 
    mutate(ral_id = gsub("DGRP", "line", ral_id))
  dim(dt)
  
  
  
  
  
  #re organize
  colnames(dt)[3] <- "avg"
  #
  #adjust the line means by inversion and wolbachia status
  #create a model that sees pheno values as a function of these covariates
  
  model = lm(avg ~ Inversion_2L_t_NA + Inversion_2R_NS_NA + Inversion_3R_K_NA + Inversion_3R_Mo_NA + Inversion_3R_P_NA + WolbachiaStatus_NA, data = dt)
  residuals = resid(model)
  dt$residuals = residuals
  #check that all dgrp lines used are in the dgrp
  succeed = dt$ral_id %in% colnames(grm)#soem number issues: 31 in grm/gds is 031
  #try making to numberic
  dt$ral_id = gsub( "line_","", dt$ral_id)
  dt$ral_id2 = as.numeric(dt$ral_id)#numeric removes excess zeroes
  dt$ral_id2 = paste0("line_", dt$ral_id2)
  #check against null grm
  
  #GMMAT creates a null model that uses wolbachia status as fixed effect on phenotype, and the identity matrix as a (nonexistent) random effect. family refers to the method used to calculate the model
  #fit a GLMM to the data
  modelqtl <- glmmkin(fixed = residuals ~ 1, data = dt, kins = grm, id = "ral_id2",family = gaussian(link = "identity"))
  
  
  
  #run GMMAT
  #create  path for output gmmat file to be saved to- using the phenotype variable allows it to be saved to the corresponding directory
  outputs = paste(INPUT , 
                  "permOctfullGRMGWAS/", 
                  f,"-",pheno.id, ".txt",
                  sep = "")
  
  print(str(modelqtl))
  print(ingds)
  print(outputs)
  
  
  #the glmm.score function scores the impact of each snp agains the null model. minor allele frequency and missing data cutoff can be specified to reduce unwanted computation. compututation can also be sped up by specifying the number of snps calculated per batch, and how many cores to use
  glmm.score(modelqtl, infile = ingds, outfile = outputs, MAF.range = c(0.05,0.95), miss.cutoff = 0.15, ncores = numCores ,
             nperbatch = 100, verbose = T)
  
  
  print ("done")
  # #replace file with rds file to save space
  # x = fread(outputs)
  # saveRDS(x, paste(INPUT , 
  #                  "OctcovarGRMGWAS/", 
  #                  phenolist[id],
  #                  sep = ""))
  # file.remove(outputs)#remove txt file
  # print("finished")
  
}

