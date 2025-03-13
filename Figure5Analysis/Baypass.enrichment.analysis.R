#all enrichments processing
library(data.table)
library(foreach)
library(doMC)
registerDoMC(5)
library(tidyverse)
setwd("/scratch/bal7cg/fall2024objects/")


b.data =readRDS("baypass.summary.oct")
keeptraits = fread("keeptraits.csv", header = T)
sig.phenos = readRDS("sigphenos")

##########################################################
#see about binding the two files with or as a trait
#make two combined files- one for chr/inv, one without
#chr inverted file.

b.cleaned = b.data %>% 
  select(-c(TTb,TFb,FTb,FFb,TTx,TFx,FTx,FFx, or.b.y, or.x.y)) %>% 
  pivot_longer(c(or.b.x, or.x.x), names_to = "type", values_to = "or"
  )

##########################################
#find quantile range for permutatations

all.data = b.cleaned %>% 
  # pivot_longer(c(or.e, or.c), names_to = "sample.set", values_to = "or") %>%  
  mutate(method = case_when(method == "OctcovarGRMGWAS" ~ "Factored.out",
                            method == "OctFullGRMGWAS" ~ "FullGrm",
                            method == "OctLDGRMGWAS" ~ "LD",
                            method == "OctLocoGRMGWAS" ~ "Loco",
                            method == "OctNoGRMGWAS" ~ "None")) %>% 
  mutate(type = case_when(type == "or.b.x" ~ "bayes.enrichment",
                           type == "or.x.x" ~ "XtX.enrichment",
                           T ~ "annotation"))
testdata = all.data %>%
  #filter(perm.number != 0) %>% 
  group_by(chr,inv.st, method,  phenotype,perm.st, threshold, type) %>%
  summarise( lowerlimitor = quantile(or, 0.05, na.rm = T),
             upperlimitor = quantile(or, 0.95, na.rm = T)) %>%
  as.data.table(.)

#####################
#rank observed data
obs.dt = testdata[perm.st == "observed"]
perm.dt = testdata[perm.st == "permutation"]

#create an object of 
obs.dt = obs.dt %>% 
  group_by(chr,inv.st, method,perm.st, threshold, type) %>% 
  mutate(rank = rank(lowerlimitor)) %>% 
  as.data.table(.)
#try merging in perm data based on pheno
mergedata = perm.dt %>% 
  select( inv.st, chr, phenotype, threshold, type, method , lowerlimitor, upperlimitor) %>% 
  merge(obs.dt, ., by = c("inv.st", "chr", "phenotype", "threshold", "type", "method")) %>% 
  mutate(sig = case_when(lowerlimitor.x > upperlimitor.y ~ "Sig",
                         
                         T ~ "Non-Sig")) 
saveRDS(mergedata, "east.coast.chr.ranked")
