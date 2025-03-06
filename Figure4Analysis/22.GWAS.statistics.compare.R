#library(GMMAT)
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gmodels)
registerDoParallel(4)

#check out results
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Sep_2023_objects//")
gifdata = readRDS("gwas.bind")
#key take away: we loaded in each of the loco files, but because they are chromosome specific, they would not work with the non-matching chromsoome sort. otherwise they should be ok, just filter out with na omit
# gifdata = gifdata[grepl("_2L", gifdata$pheno) == T ,]
# gifdata = gifdata[grepl("_2R", gifdata$pheno) == F ,]
# gifdata = gifdata[grepl("_3L", gifdata$pheno) == F ,]
# gifdata = gifdata[grepl("_3R", gifdata$pheno) == F ,]
# gifdata = gifdata[grepl("_X", gifdata$pheno) == F ,]
dim(is.na(gifdata))
gifclean = na.omit(gifdata)
#remove chromosome names from loco phenos
gifclean$pheno = gsub("_2L", "", gifclean$pheno)
gifclean$pheno = gsub("_2R", "", gifclean$pheno)
gifclean$pheno = gsub("_3L", "", gifclean$pheno)
gifclean$pheno = gsub("_3R", "", gifclean$pheno)
gifclean$pheno = gsub("_X", "", gifclean$pheno)

#hype!
####analysis###
#we want to see how the metrics change across inversion, chrom, and grm. it would be cool to eventually add pheno metadata
#scale sig snps
gifclean$sigratio = gifclean$Sig.number / gifclean$Totalsnpnumber

#load in perm data objects
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/May_2024_objects//")
permfull = get(load("permgather,1"))
permfull$grm = "OctFullGRMGWAS"
permld = get(load("permgather,2"))
permld$grm = "OctLDGRMGWAS"
permcovar=get(load("permgather.covar"))
permcovar$grm = "OctcovarGRMGWAS"
# permcovar2 = get(load("permcovar"))
# permcovar2$grm = "OctcovarGRMGWAS"
permno=get(load("permgather,1-no"))
permno$grm = "OctNoGRMGWAS"
permloco=get(load("permgatherloco2"))
permloco$grm = "OctLocoGRMGWAS"
permdata = rbind(permfull, permld, permno, permloco,permcovar)
permdata$sigratio = permdata$Sig.number / permdata$Totalsnpnumber
#what if we try to find the mean gif per permutation? 
# permdata[,permnumber := tstrsplit(pheno, "-")[1]]
# permdata = permdata %>% 
#   mutate(perm.number = )
#   group_by(chromosome, inversion, grm, )

gifclean$perm.st = "observed"
#gifclean$permnumber = NA
permdata = na.omit(permdata)
permdata$perm.st = "permutations"
colnames(gifclean)
colnames(permdata)
#before we bind, we need to remove perm numbers and file stems from perm data
permdata[,pheno := tstrsplit(permdata$pheno, "-")[2]]
permdata$pheno = gsub(".txt .RDS", "", permdata$pheno)
permdata$pheno = gsub(".txt", "", permdata$pheno)#ok, we see that some perms did not fully finish their gwas(still in gmmats tmp form)
permdata2 = permdata[grepl( "tmp", permdata$pheno) == F]
#check how many phenos succesffuly gathered per method
perm.count = permdata2 %>% 
  group_by(grm) %>% 
  summarise(pheno) %>% 
  distinct(.) %>% 
  summarise(N = n())
#it seems loco and no grm messed up - no phenotype is recorded
#covar is still only gathering half of the phenotypes
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/July_2024_objects/")
fulldata = rbind(gifclean, permdata2)
fulldata = fulldata %>% 
  mutate(grm = case_when(grm == "OctcovarGRMGWAS" ~ "Factored.out" ,
                         grm == "OctFullGRMGWAS" ~ "FullGrm",
                         grm == "OctLDGRMGWAS" ~ "LD",
                         grm == "OctLocoGRMGWAS" ~ "Loco",
                         grm == "OctNoGRMGWAS" ~ "None"))
x =fulldata[grm == "None"]
saveRDS(fulldata, "gwas.bind.data.aug")
#lets make plots showing distriution of gif and #sig phenos for each group /real

dis.data = fulldata %>% 
  #filter(grm == "full") %>% 
  group_by(chromosome, inverted, pheno, perm.st,grm) %>%      summarise(avg = mean(sigratio),
                                                                        lowerbound = quantile(sigratio, 0.025),
                                                                        upperbound = quantile(sigratio, 0.975)) %>% 
  as.data.table(.)
saveRDS(dis.data, "sig.grouped.data.aug")
#now somehow make the graph- start with the observed.
dis.data$pheno = gsub(".RDS", "", dis.data$pheno)
obs.dt = dis.data[perm.st == "observed"]

perm.dt = dis.data[perm.st == "permutations"]
#we want to make a seperate column that ranks the phenos based on their gif observed, then match in their perm values
obs.dt = obs.dt %>% 
  group_by(inverted, chromosome,grm) %>% 
  mutate(rank = rank(avg)) %>% 
  as.data.table(.)
#try merging in perm data based on pheno
mergedata = perm.dt %>% 
  select(chromosome, inverted, pheno, avg, lowerbound, upperbound, grm) %>% 
  merge(obs.dt, ., by = c("inverted", "chromosome", "pheno", "grm")) %>% 
  mutate(sig = case_when(avg.x > upperbound.y ~ "Sig",
                         avg.x < lowerbound.y ~ "Sig",
                         T ~ "Non-Sig")) %>% 
  as.data.table()
saveRDS(mergedata, "snpperm.obs.compare")
# #use this source script to try reordering
# source("/Users/supad/OneDrive/Documents/Bergland Research/R_scripts/May_2024.scripts/reorder_within.R")
# x = dis.data %>% 
#   group_by(pheno) %>% 
#   summarise(mean.avg = mean(avg)) %>% 
#   arrange(mean.avg) 
full.dt = mergedata[grm == "None"]


ggplot() +
  geom_point(data = full.dt, aes(x = rank, y = avg.x, color = sig))+
  geom_errorbar(data = full.dt, aes(x=rank,
                                  y=avg.y,
                                  ymin=lowerbound.y,
                                  ymax=upperbound.y),
                             width = 0.1, position=position_dodge(width = 0.5), show.legend = F, alpha = 0.05) +
  facet_grid(inverted~ chromosome)+
  xlab("phenotypes")+
  ylab("GIF")+
  #ggtitle("FullGRM Sig SNP Distributions")+
  theme_bw()
             
             geom_point(position=position_dodge(width = 0.5), show.legend = T)