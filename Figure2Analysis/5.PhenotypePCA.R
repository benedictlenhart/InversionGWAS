#run PCA analysis on phenotypes associated with inversions, and get loading info for graphing

setwd( wd )
### libraries
library(data.table)
library(ggplot2)
library(SeqArray)
library("Hmisc")
library(FactoMineR)
library(factoextra)
library(missMDA)
library(tidyverse)
library(patchwork)
library(foreach)
### load
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Sep_2024_objects/")
#load in phenotype data
pheno <- readRDS("sepphenolongform")
dim(pheno)


setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2025_objects/")
#load in identifies of the inversion associated traits
sig.phenos = readRDS("C:/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Feb_2025_objects/topmost.phenos")
#############################################
###compare the two sets of sig phenos#########
############################################3

inversions = unique(sig.phenos$inversion)
#ref to the column we want from the inv file
loading.list = list()
pca.list = list()
#make a foreach loop that goes through each inversion

c.out = foreach(f = c(1:5)) %do%{
  # f = 1
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
  
  #another format of inversion genotype info
  ### inversion
  inv <- fread("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Sep_2024_objects//inversion.csv", header = T)
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
  imputed = imputePCA(pw[,-c("DGRP"), with=F])
  imputed.table = as.data.frame(imputed[1])
  if (f == 5 ){
    pc <- PCA(imputed, scale.unit=T, graph=F)
  }else {
    pc <- PCA(imputed.table, scale.unit=T, graph=F)
  }
  
  #pc = PCA(pheno, scale.unit = T, graph = F)
  
  #viz <- fviz_pca_var(pc, select.var = list(name = NULL, cos2 = NULL, contrib = NULL), repel=F) + theme(text=element_text(size=1))
  
  
  pc_loading <- as.data.table(pc$var$coord)
  pc_loading[,pheno:=rownames(pc$var$coord)]
  step=.1
  ### top left quadrant
  pc_loading[Dim.1<0 & Dim.2>0, x:=-.4]
  pc_loading[Dim.1<0 & Dim.2>0, 
             y:=1.5 - seq(from=0, length.out=dim(pc_loading[Dim.1<0 & Dim.2>0])[1], by=step)]
  
  ### bottom left quadrant
  pc_loading[Dim.1<0 & Dim.2<0, x:=-.25]
  pc_loading[Dim.1<0 & Dim.2<0, 
             y:=-1 - seq(from=0, length.out=dim(pc_loading[Dim.1<0 & Dim.2<0])[1], by=step)]
  
  ### bottom right quadrant
  pc_loading[Dim.1>0 & Dim.2<0, x:=.25]
  pc_loading[Dim.1>0 & Dim.2<0, 
             y:=-.5 - seq(from=0, length.out=dim(pc_loading[Dim.1>0 & Dim.2<0])[1], by=step)]
  
  ### top right quadrant lower on Dim1 (TOP CLUSTER)
  pc_loading[Dim.1>0 & Dim.1 <=0.25 & Dim.2>0, 
             x:=.05]
  pc_loading[Dim.1>0 & Dim.1 <=0.25 & Dim.2>0, 
             y:=1.5 - seq(from=0, length.out=dim(pc_loading[Dim.1>0 & Dim.1 <=0.25 & Dim.2>0])[1], by=step)]
  
  ### top right quadrant upper on Dim1
  pc_loading[Dim.1>.25 & Dim.2>0, x:=.8]
  pc_loading[Dim.1>.25 & Dim.2>0, 
             y:=.65 - seq(from=0, length.out=dim(pc_loading[Dim.1>.25 & Dim.2>0])[1], by=step)]
  
  ### indiviuals PCA data table
  pcr <- cbind(pw[,c("DGRP"), with=F], as.data.table(pc$ind$coord))
  #filter 
  pcr <- merge(pcr, inv[V2 != "INV/ST"], by="DGRP") #remove heterozygotes for now
  
  # summary(lm(Dim.1~gt, pcr))
  # summary(lm(Dim.2~gt, pcr))
  # summary(lm(Dim.3~gt, pcr))
  
  summary(lm(Dim.1~V2, pcr))
  summary(lm(Dim.2~V2, pcr))
  summary(lm(Dim.3~V2, pcr))
  
  
  
  
  pcr.ag <- pcr[!is.na(V2),list(mu1=mean(Dim.1), se1=sd(Dim.1)/sqrt(length(Dim.1)),
                                mu2=mean(Dim.2), se2=sd(Dim.2)/sqrt(length(Dim.2)),
                                mu3=mean(Dim.3), se3=sd(Dim.3)/sqrt(length(Dim.3))),
                list(V2)]
  pcr.ag[V2== "INV", gt.name:="Inverted"]
  pcr.ag[V2== "ST", gt.name:="Standard"]
  
  ### loadings plot
  pcl <- pc$var$coord
  pclw <- as.data.table(melt(pcl))
  pclw[,dim:=as.numeric(tstrsplit(Var2, "\\.")[[2]])]
  pcr.ag$inv = inv.of
  pc_loading$inv = inv.of
  loading.list= list(pcr.ag, pc_loading)
  loading.list
}
#okay so it's a list of five smaller lists- how to get them out?
loadings = list( c.out[[1]][[1]], c.out[[2]][[1]], c.out[[3]][[1]],c.out[[4]][[1]],c.out[[5]][[1]])
loadings = rbindlist(loadings)
pcas = list( c.out[[1]][[2]], c.out[[2]][[2]], c.out[[3]][[2]],c.out[[4]][[2]],c.out[[5]][[2]])
pcas = rbindlist(pcas, fill = T)

saveRDS(loadings, "pca_loading.data2")
saveRDS(pcas, "pca_dim.data2")
