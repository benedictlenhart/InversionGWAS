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
#id = 1
#id = 1#use this when testing to see if script functions correctly


# saveRDS(task.id, paste0(INPUT, "newphenotaskidchart"))
INPUT = "/project/berglandlab/Adam/fall2023objects/"#path of working directory
ingds = paste0(INPUT,"dgrp.gds")
#load in reftable
reftable = readRDS(paste0(INPUT,"dgrpreftable"))
#load grm
grm = readRDS(paste0(INPUT, "GRM.Aug"))
#load in our phenotype file 
#phenotable = readRDS(paste0(INPUT, "sepphenolongform"))
phenotable = readRDS(paste0(INPUT, "pca.long.imputed4"))
#remake pca data to fit older format
phenotable = phenotable %>% 
  mutate(fullpheno = paste(Inversions, Dim, sep = "-")) %>% 
  rename("avg" = "loadings") %>% 
  as.data.table(.)
#phenotable = readRDS(paste0(INPUT, "sepphenowideform"))
#filter down to phenotype of interet
phenolist = unique(phenotable$fullpheno)


dt = phenotable[fullpheno == phenolist[id]]
#merge in the values for reftable
dt = merge(dt, reftable, by = "DGRP" )
dt = na.omit(dt)#remove na's 
#fix ral ids
dt = dt %>% 
  dplyr::rename("ral_id" = DGRP) %>% 
  mutate(ral_id = gsub("DGRP", "line", ral_id))
dim(dt)





#re organize

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
                "factored_out.pca.gwas2/", 
                phenolist[id], ".txt",
                sep = "")


#the glmm.score function scores the impact of each snp agains the null model. minor allele frequency and missing data cutoff can be specified to reduce unwanted computation. compututation can also be sped up by specifying the number of snps calculated per batch, and how many cores to use
glmm.score(modelqtl, infile = ingds, outfile = outputs, MAF.range = c(0.05,0.95), miss.cutoff = 0.15, 
           nperbatch = 400, ncores = 5)


print ("done")
#replace file with rds file to save space
x = fread(outputs)
saveRDS(x, paste(outputs,
                 ".RDS"))
file.remove(outputs)#remove txt file
print("finished")

