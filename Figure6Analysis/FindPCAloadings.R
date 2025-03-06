#system("cp /project/berglandlab/Yang_Adam/reference_files/dgrp.gds .")


### libraries
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

# new.sig.phenos = readRDS("sigphenos.frommodeling"

inversions = unique(sig.phenos$inversion)
#ref to the column we want from the inv file
loading.list = list()
pca.list = list()
#make a foreach loop that goes through each inversion

c.out = foreach(f = c(1:5)) %do%{
 #f = 3
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
    # mutate(Inversions = case_when(Inversions == "In(2L)t" ~ "2L_t",
    #                               Inversions == "In(2R)NS" ~ "2R_NS",
    #                               Inversions == "In(3R)P" ~ "3R_P",
    #                               Inversions == "In(3R)Mo" ~ "3R_Mo",
    #                               Inversions == "In(3R)K" ~ "3R_K",
    #                               T ~ "other"))%>% 
    filter(Inversions == inv.of) %>% 
    as.data.table(.)
  
  ### principal components
  #foreach(unique(pl$pos))%do%{
  #now
  #pw <-dcast(pl[pos==5192177][variable%in%target$pheno], ral_id+gt~variable, value.var="x_norm")
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
  ### indiviuals PCA data table
  pcr <- cbind(pw[,c("DGRP"), with=F], as.data.table(pc$ind$coord))
  #filter 
  pcr <- merge(pcr, inv[V2 != "INV/ST"], by="DGRP") #remove heterozygotes for now
  #same the pcas
  pcr$Inversions = inv.of
  pcr
}
pcas = rbindlist(c.out, fill = T)
#make data longform
pca.long = pcas %>%
  pivot_longer(!c(DGRP, Inversions, V2), names_to = "Dim", values_to = "loadings") %>%
  rename("Inv.st"="V2")
# pca.long = pcas %>%
pca.dt = pca.long
#   pivot_longer(!c(DGRP,Inversions), names_to = "Dim", values_to = "loadings") 
  saveRDS(pca.long, "pca.long.imputed4")
pca.dt =  readRDS("pca.long.imputed4")
##############################################################
###find which pcas are significantly impacted by inversion###
#############################################################
sigs = pca.dt %>% 
  group_by(Dim, Inversions) %>% 
  summarize(pvalue = t.test(loadings ~ Inv.st, data = .)$p.value) %>% 
  mutate(is.sig = ifelse(pvalue < 0.05, "Sig", "Not-Sig")) %>% 
  arrange(Inversions)
  pca.dt = as.data.table(pca.long)
  ref.table = expand.grid(unique(pca.dt$Dim),unique(pca.dt$Inversions))
sigs
out = foreach(f = c(1:dim(ref.table)[1]), .errorhandling = "remove") %do% {
# f =1
  ref.info = ref.table[f,]
  pca.dt = as.data.table(pca.dt)
  dt = pca.dt[Dim== unlist(ref.info[1,1])][Inversions == unlist(ref.info[1,2])]
  # model.add = lm(value ~ def.id + inv.st , data = dt)
  # model.mult = lm(value ~ def.id * inv.st , data = dt)
  # afit = anova(model.add, model.mult)
  I.S = t.test(dt[Inv.st == "INV"]$loadings,dt[Inv.st == "ST"]$loadings)
  
  dt = data.frame(
    Dim = unlist(ref.info[1,1]),
    Inversions = unlist(ref.info[1,2]),
  
    df_I.S = I.S$parameter,
    t_I.S = I.S$statistic,
    p_I.S = I.S$p.value

  )
  dt
}
stats.out = rbindlist(out)
#re
saveRDS(stats.out, "/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Feb_2025_objects/stats.out")
tall = stats.out %>% 
  pivot_longer(cols = c(p_I.S, p_H.S,p_H.I ), names_to = "comparisons")
sig = tall %>% 
  filter(value < 0.05) %>% 
  as.data.table(.)
sig

  