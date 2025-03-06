#fig 1 linear modeling graphing

library(lme4)
library(data.table)
library(tidyverse)
library(foreach) 

setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2025_objects/")

fulldata = readRDS("same.freq.obj2")

keeptraits = readRDS("keeptraits.fixed")




#filter out the inversion name from inversions
fulldata$inversion = gsub("Inversion_", "", fulldata$inversion)
fulldata$inversion = gsub("_NA", "", fulldata$inversion)
#filter to our final phenotypes and add meta info

fulldata = fulldata %>%

  merge(keeptraits, ., by = "phenotype") %>%
  mutate(meta.info = case_when(meta.info == "p"~ "Physiology",
                               meta.info == "sr"~ "Stress_Resistance",
                               meta.info == "lh"~ "Life_History",
                               meta.info == "m"~ "Morphology",
                               meta.info == "b"~ "Behavior")) %>%
  as.data.table(.)
sigphenos = fulldata %>%
  filter(P.value < 0.05,
         snp.pos == "inversion")

#find which results are significant
sigresults = fulldata %>%
  mutate(snp.pos = case_when(snp.pos == "inversion" ~ "0",
                             T ~ snp.pos),
         snp.pos = as.numeric(snp.pos))%>%
  filter(snp.pos < 101) %>%
  # filter(P.value < 0.05) %>%
  mutate(perm.st = case_when(snp.pos == 0 ~ "observed",
                             T ~ "Same.Freq.Snps")) %>%
  group_by(inversion, snp.pos,perm.st, meta.info, meta.count) %>%
  #filter(P.value < 0.05) %>%
  summarise(number = sum(P.value < 0.05, na.rm = T)) %>%
  as.data.table(.)

sigresults= sigresults %>% 
  mutate(inversion = case_when(inversion == "2L_t"~"In(2L)t",
                               inversion == "2R_NS"~"In(2R)NS",
                               inversion == "3R_K"~"In(3R)K",
                               inversion == "3R_Mo"~"In(3R)Mo",
                               inversion == "3R_P"~"In(3R)P"))

##############
##graphing###
#############

#create a plot that shows perms as grey, real as red
g1 =ggplot() +
  geom_boxplot(width = 0.5, data = sigresults[perm.st == "Same.Freq.Snps"], 
               aes(x = inversion, y = number, color = meta.info), outlier.colour = "NA", show.legend = F, position =position_dodge(width= 0.5)) +
  # geom_jitter(data =sigresults[perm.st == "Same.Freq.Snps"], 
  #             aes(x = inversion, y = number, color = meta.info), alpha = 0.08, width = .2, show.legend = F) +
  ylim(0,23)+
  geom_point(data = sigresults[perm.st == "observed"],
             shape = 23, 
             size = 5,
             position = position_dodge(width= 0.5),
             aes(color = meta.info,
                 group = meta.info,
                 fill = meta.info,
                 x = inversion,
                 y = number),
             show.legend = T
  ) +
  #facet_grid(cols = vars(inv), scales = "free")+
  theme_bw(base_size = 15) + 
  #theme(text = element_text(size = 13))+
  #facet_grid(method~.) +
  xlab("Inversions") +
  ylab( "# of Sig. Phenotypes (< 0.05) ")
g1

#Next, graph distribution of r squared
#will need to rank each phenotype by highest r square.

mean.data = fulldata %>%
  mutate(snp.pos = case_when(snp.pos == "inversion" ~ "0",
                             T ~ snp.pos),
         snp.pos = as.numeric(snp.pos))%>%
  filter(snp.pos < 101) %>%
  na.omit(.) %>%
  mutate(perm.st = case_when(snp.pos == "0" ~ "observed",
                             T ~ "Same.Freq.Snps")) %>%
  group_by(phenotype, snp.pos, inversion, perm.st) %>%
  summarise(mean.r = mean(rsquared, na.rm = T))
# saveRDS(mean.data , "mean.data.mod")
# mean.data = readRDS("mean.data.mod")

dis.data = mean.data%>% 
  
  group_by( inversion, phenotype, perm.st) %>%      
  summarise(avg = mean(mean.r, na.rm = T),
            lowerbound = quantile(mean.r, 0.025),
            upperbound = quantile(mean.r, 0.975)) %>% 
  mutate(inversion = case_when(inversion == "2L_t"~"In(2L)t",
                               inversion == "2R_NS"~"In(2R)NS",
                               inversion == "3R_K"~"In(3R)K",
                               inversion == "3R_Mo"~"In(3R)Mo",
                               inversion == "3R_P"~"In(3R)P")) %>% 
  as.data.table(.)
obs.dt = dis.data[perm.st == "observed"]
perm.dt = dis.data[perm.st == "Same.Freq.Snps"]
sigphenos = sigphenos %>% 
  #rename("Classification"="meta.info") %>% 
  mutate(inversion = case_when(inversion == "2L_t"~"In(2L)t",
                               inversion == "2R_NS"~"In(2R)NS",
                               inversion == "3R_K"~"In(3R)K",
                               inversion == "3R_Mo"~"In(3R)Mo",
                               inversion == "3R_P"~"In(3R)P")) %>% 
  select(inversion, phenotype)
obs.dt = obs.dt %>% 
  merge(., sigphenos, by = c("inversion","phenotype")) %>% #only graph our linear sig phenos
  group_by(inversion) %>% 
  mutate(rank = rank(avg)) %>% 
  as.data.table(.)
#try merging in perm data based on pheno
mergedata = perm.dt %>% 
  select( inversion, phenotype, avg, lowerbound, upperbound) %>% 
  merge(obs.dt, ., by = c("inversion", "phenotype")) %>% 
  mutate(sig = case_when(avg.x > upperbound.y ~ "Sig",
                         
                         T ~ "Non-Sig")) %>% 
  # filter(avg.x <= 0.25) %>% 
  as.data.table()
#make a data table that counts number of sig results and can be added in as labels
count = mergedata %>% 
  group_by(inversion) %>% 
  filter(sig == "Sig") %>%
  #filter(rank > 250) %>% 
  summarise(count = n()) %>% 
  mutate(rank = 10,
         avg.x = .3)

colors = c("black","cyan")
g2 =ggplot() +
  geom_point(data = mergedata, aes(x = rank, y = avg.x, color = sig),show.legend = T)+
  geom_errorbar(data = mergedata, aes(x=rank,
                                      y=avg.y,
                                      ymin=lowerbound.y,
                                      ymax=upperbound.y),
                width = 0.1, position=position_dodge(width = 0.5), show.legend = F, alpha = 0.2) +
  facet_grid(.~inversion, scales = "free")+
  xlab("Phenotypes")+
  ylab(bquote(R^2))+
  scale_x_continuous(n.breaks = 3)+
  scale_color_manual(values = colors)+
  ylim(0, .75)+
  #ggtitle("FullGRM Sig SNP Distributions")+
  theme_bw(base_size = 15) +
  geom_text(data = count,aes(x = rank, y = avg.x, label = count)) 
# + 
#   theme(text = element_text(size = 13))

g2

#find super sig traits- those that are significant from linear modeling and surpass the same.freq model
#these are the traits we are determing as inversion-linked groups
sigphenos$merge = paste(sigphenos$phenotype, sigphenos$inversion, sep = "-")

mergedata$merge = paste(mergedata$phenotype, mergedata$inversion, sep = "-")
top.sig = mergedata %>% 
  filter(sig == "Sig") %>%
  merge(., sigphenos, by = c("phenotype", "inversion")) %>% 
  select(c(phenotype, inversion))
saveRDS(top.sig, "C:/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Feb_2025_objects/topmost.phenos")
####
#save graphs
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2025_objects/images/")
library(patchwork)
ggsave((( g1 + g2 ) + plot_layout(guides = "collect") & theme(legend.position = 'bottom'))+
         plot_annotation(tag_levels = 'A'), file = "gwasfig1.5.pdf", height = 5, width = 13) 

#################################################
###add significance data to our meta data sheet##
#################################################
sigphenos = fulldata %>%
  filter(P.value < 0.05,
         snp.pos == "inversion")
meta = fread("Supp.meta.info.csv")

#create column for meta that is linear model signficance for each inversion
sigphenos = sigphenos %>% 
  mutate(inversion = case_when(inversion == "2L_t"~"In(2L)t",
                               inversion == "2R_NS"~"In(2R)NS",
                               inversion == "3R_K"~"In(3R)K",
                               inversion == "3R_Mo"~"In(3R)Mo",
                               inversion == "3R_P"~"In(3R)P"))
inv = unique(sigphenos$inversion)
#now create a foreach loop that filters delta and chr.ref by chr, then joins
chr.out = foreach(f = inv) %do% {
  #f = "In(2L)t"
  chr.phenos = sigphenos%>% 
    filter(inversion == f) 
  chr.phenos = unlist(chr.phenos$phenotype)
  meta2 <- function(df, f) {
    newcol = paste("linear_model:", f, sep = "_")
    df[[newcol]] <- meta$fullpheno %in% chr.phenos
    df
  }
  meta = meta2(meta, f)
  
}
#now repeat with r2 to show our sigphenos
sig.phenos = readRDS("C:/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Feb_2025_objects/topmost.phenos")
chr.out = foreach(f = unique(sig.phenos$inversion)) %do% {
  #f = "In(2L)t"
  chr.phenos = sigphenos%>% 
    filter(inversion == f) 
  chr.phenos = unlist(chr.phenos$phenotype)
  meta2 <- function(df, f) {
    newcol = paste("Significantly_associated:", f, sep = "_")
    df[[newcol]] <- meta$fullpheno %in% chr.phenos
    df
  }
  meta = meta2(meta, f)
  
}
write.csv(meta, "Supp.dataV2.csv")