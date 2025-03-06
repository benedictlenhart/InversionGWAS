library(data.table)
library(ggplot2)
library(SeqArray)
library("Hmisc")
library(FactoMineR)
library(missMDA)
library(factoextra)
library(tidyverse)
library(patchwork)
library(foreach)
### load
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Sep_2024_objects/")
pheno <- readRDS("sepphenolongform")

setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2025_objects/")

sig.phenos = readRDS("C:/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Feb_2025_objects/topmost.phenos")


inversions = unique(sig.phenos$inversion)
#ref to the column we want from the inv file
loading.list = list()
pca.list = list()
#make a foreach loop that goes through each inversion

c.out = foreach(f = c(1:5)) %do%{
 # f = 3
  inv.of = inversions[f]
  #we want to filter to only phenos of interest (say in2lt, loco method)
  #merge pheno and inversion info
  #pheno = merge(pheno, dgrpref, by = "DGRP")
  interest = sig.phenos %>% 
    filter(inversion == inv.of) %>% 
    select(phenotype) %>% 
    unlist(.)
  
  #narrow phenodata to only phenotyps of interst
  narrow = pheno %>% 
    filter(fullpheno %in% interest) %>% 
    select(!phenotype.sex)
  #now narrow ref to only inversion of interest
  
  
  ### inversion
  inv <- fread("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Sep_2024_objects/inversion.csv", header = T)
  inv$DGRP = inv$`DGRP Line`
  #melt and filter the data
  inv = inv %>% 
    pivot_longer(
      !DGRP,
      names_to = "Inversions",
      values_to = "V2"
    ) %>% 
  
    filter(Inversions == inv.of) %>% 
    as.data.table(.)
  
  ### principal components
  
  #for this to work need narrow to have dgrp ides, and all other phenos
  pw = dcast(na.omit(narrow),DGRP ~ fullpheno, value.var = "pheno.mean")
  #use pca based imputation to fill in missing data
  imputed = imputePCA(pw[,-c("DGRP"), with=F])
  imputed.table = as.data.frame(imputed[1])
  #w <-dcast(pl, ral_id+gt~variable, value.var="x_norm")
  
  if (f == 5 ){
    pc <- PCA(imputed, scale.unit=T, graph=F)
  }else {
    pc <- PCA(imputed.table, scale.unit=T, graph=F)
  }
  
  #get the eigenvalues for each pc.
  out = data.frame(
    eigen = as.numeric((pc$eig[,2]))
  )
  out$PC = c(1:dim(out)[1])
  out$Inversion = inv.of
  out
}
pcas = rbindlist(c.out, fill = T)
saveRDS(pcas,"scree.var")
