---
  title: "gwas.fig4"
output: html_document
date: "2024-08-07"
---
#15 minutes into chapte r6
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

libraries and data
```{r}
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gmodels)
registerDoParallel(4)

#check out results
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/July_2024_objects/")
gifdata = readRDS("perm.obs.compare")
snpdata = readRDS("snpperm.obs.compare")
loppdata = readRDS("lowp.obs.compare")
#testanno = readRDS("enrichment_OctLocoGRMGWAS-AbdominaPigmentT5_F")
#load in keeppheno names
#check number of traits passing 

setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2025_objects/")
keeptraits = readRDS("keeptraits.fixed")
sigphenos = readRDS("C:/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Feb_2025_objects/topmost.phenos")




# Alan suggests a slight modication where we get "p", the proportion of phenotypes that beat perms, then use that generate binoial confidence interavals with formula mean +/-  1.96 * sqrt(p(1-p)/n)


#we need to get the proportion instead of the count, then use that generate upper and lower intervals
propgif1 = gifdata %>% 
  rename("phenotype"= "pheno") %>% 
  merge(keeptraits, ., by = "phenotype" ) %>% 
  group_by( inverted, chromosome, grm) %>% 
  mutate(totalcount = n())
propgif2 = propgif1 %>% 
  filter(sig == "Sig") %>% 
  filter(avg.x > upperbound.y) %>% #we only want right tail
  group_by( inverted, chromosome, grm, totalcount) %>% 
  summarize(pass.count = n()) %>% 
  mutate(p = pass.count / totalcount,
         upperci = p + 1.96 * sqrt((p *(1-p)) / totalcount),
         lowerci = p - 1.96 * sqrt((p *(1-p)) / totalcount),
         trait = "GIF")

propsnp1 = snpdata %>% 
  rename("phenotype"= "pheno") %>% 
  merge(keeptraits, ., by = "phenotype" ) %>% 
  group_by(inverted,chromosome, grm) %>% 
  mutate(totalcount = n())
propsnp2 = propsnp1 %>% 
  filter(sig == "Sig") %>% 
  filter(avg.x > upperbound.y) %>% #we only want right tail
  group_by(inverted, chromosome, grm, totalcount) %>% 
  summarize(pass.count = n()) %>% 
  mutate(p = pass.count / totalcount,
         upperci = p + 1.96 * sqrt((p *(1-p)) / totalcount),
         lowerci = p - 1.96 * sqrt((p *(1-p)) / totalcount),
         trait = "SNpcount")

propp1 = loppdata %>% 
  rename("phenotype"= "pheno") %>% 
  merge(keeptraits, ., by = "phenotype" ) %>%  
  filter(chromosome != 4) %>% 
  group_by(inverted, chromosome, grm) %>% 
  mutate(totalcount = n())
propp2 = propp1 %>% 
  filter(sig == "Sig")%>% 
  filter(avg.x < lowerbound.y) %>% #we only want right tail
  group_by(inverted, chromosome, grm, totalcount) %>% 
  summarize(pass.count = n()) %>% 
  mutate(p = pass.count / totalcount,
         upperci = p + 1.96 * sqrt((p *(1-p)) / totalcount),
         lowerci = p - 1.96 * sqrt((p *(1-p)) / totalcount),
         trait = "LowestPvalue")
fullprop2 = rbind(propgif2, propsnp2,propp2) 
#simplify names
fullprop2 = fullprop2 %>% 
  mutate(inverted = gsub("inverted", "Inverted", inverted),
         inverted = gsub("non", "Non", inverted)) %>% 
  mutate(grm = case_when(grm == "OctcovarGRMGWAS" ~ "Factored-out" ,
                         grm == "Factored.out" ~ "Factored-out" ,
                         grm == "OctFullGRMGWAS" ~ "Full",
                         grm == "FullGrm" ~ "Full",
                         grm == "OctLDGRMGWAS" ~ "LD",
                         grm == "LD" ~ "LD",
                         grm == "OctLocoGRMGWAS" ~ "LOCO",
                         grm == "Loco" ~ "LOCO",
                         grm == "OctNoGRMGWAS" ~ "None",
                         T ~ grm)
  ) %>% 
  mutate(trait = case_when(trait == "GIF" ~ "GIF",
                           trait == "LowestPvalue" ~ "LowestPvalue",
                           trait == "SNpcount" ~ "#ofHits")) %>% 
  filter(grm != "None",
         trait %in% c( "#ofHits", "GIF") %>% 
 
  as.data.table(.)


g1colors = c("blue","red")
# annotation_df <- data.frame(
#   chromosome = c("2L","2R", "3L", "3R","2L","2R", "3L", "3R"),
#   trait = c(rep("#ofHits",4),rep("GIF",4)),
#   end = c(rep(3.7,8)),
#   start = c(rep(0.7,8)),
#   y = c(rep(0.3, 4),0.11, 0.3, 0.3, 0.11),
#   label = c(rep("***",8))
# )


annotation_df <- data.frame(
  chromosome = c(rep(c("2L","2R", "3L", "3R"), 3)),
  #trait = c(rep("#ofHits",4)),
  end = c(rep(2,4),rep(3,4),rep(4,4)),
  start = c(rep(1,12)),
  y = c(rep(0.15, 4),rep(0.3, 4),rep(0.45, 4)),
  label = c(rep("NS./NS. ", 4), rep("NS./NS.", 4), rep("***/***", 4))
)
g1 = ggplot() +  
  ylab("Proportion exceeding permutations") +
  xlab("Methods") +
  # scale_color_manual(values = colors)+
  facet_grid(trait ~ chromosome) +
  scale_color_manual( values = g1colors) +
  geom_errorbar(data = fullprop2, aes(
    x= grm,
    y=p,
    ymin=lowerci,
    ymax=upperci,
    color = inverted),
                width = 0.1, position=position_dodge(width = 0.5), show.legend = F) +
  geom_point(data = fullprop2, aes(
    y= p,
    x=grm,
    color = inverted),position=position_dodge(width = 0.5), show.legend = T) + theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  geom_vline(xintercept = 0.05, color = "grey", linetype = "dashed")+
  geom_signif(
    data = annotation_df,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    textsize = 3, vjust = -0.02,
    tip_length = 0.01,
    manual = TRUE
  ) + ylim(0, 0.6)
g1
#look at sig phenos vs not

chr.ref = sigphenos %>% 
  select(phenotype, inversion) %>% 
  mutate(chromosome = case_when(inversion == "In(2L)t" ~ "2L", 
                                inversion == "In(2R)NS" ~ "2R", 
                                T ~ "3R"))

snp.out = foreach(f = unique(chr.ref$chromosome), .combine = "rbind") %do% {
   #f = chr.ref$chromosome[1]
  chr.phenos = chr.ref %>% 
    filter(chromosome == f) 
  chr.phenos = unlist(chr.phenos$phenotype)
  delta.filter = snpdata %>% 
    filter(chromosome == f) %>% 
    mutate(sig.pheno = pheno %in% chr.phenos)
}


propgif1 = gif.out %>% 
  rename("phenotype"= "pheno") %>% 
  merge(keeptraits, ., by = "phenotype" ) %>% 
  group_by(sig.pheno, inverted, chromosome, grm) %>% 
  mutate(totalcount = n())
propgif2 = propgif1 %>% 
  filter(sig == "Sig") %>% 
  filter(avg.x > upperbound.y) %>% #we only want right tail
  group_by(sig.pheno, inverted, chromosome, grm, totalcount) %>% 
  summarize(pass.count = n()) %>% 
  mutate(p = pass.count / totalcount,
         upperci = p + 1.96 * sqrt((p *(1-p)) / totalcount),
         lowerci = p - 1.96 * sqrt((p *(1-p)) / totalcount),
         trait = "GIF")

propsnp1 = snp.out %>% 
  rename("phenotype"= "pheno") %>% 
  merge(keeptraits, ., by = "phenotype" ) %>% 
  group_by(sig.pheno, inverted,chromosome, grm) %>% 
  mutate(totalcount = n())
propsnp2 = propsnp1 %>% 
  filter(sig == "Sig") %>% 
  filter(avg.x > upperbound.y) %>% #we only want right tail
  group_by(sig.pheno, inverted, chromosome, grm, totalcount) %>% 
  summarize(pass.count = n()) %>% 
  mutate(p = pass.count / totalcount,
         upperci = p + 1.96 * sqrt((p *(1-p)) / totalcount),
         lowerci = p - 1.96 * sqrt((p *(1-p)) / totalcount),
         trait = "SNpcount")

propp1 = p.out %>% 
  rename("phenotype"= "pheno") %>% 
  merge(keeptraits, ., by = "phenotype" ) %>%  
  filter(chromosome != 4) %>% 
  group_by(sig.pheno, inverted, chromosome, grm) %>% 
  mutate(totalcount = n())
propp2 = propp1 %>% 
  filter(sig == "Sig")%>% 
  filter(avg.x < lowerbound.y) %>% #we only want right tail
  group_by(sig.pheno, inverted, chromosome, grm, totalcount) %>% 
  summarize(pass.count = n()) %>% 
  mutate(p = pass.count / totalcount,
         upperci = p + 1.96 * sqrt((p *(1-p)) / totalcount),
         lowerci = p - 1.96 * sqrt((p *(1-p)) / totalcount),
         trait = "LowestPvalue")
fullprop2 = rbind(propgif2, propsnp2,propp2) 
fullprop2 = rbind(propsnp2)
#simplify names
fullprop2 = fullprop2 %>% 
  mutate(inverted = gsub("inverted", "Inverted", inverted),
         inverted = gsub("non", "Non", inverted)) %>% 
  mutate(grm = case_when(grm == "OctcovarGRMGWAS" ~ "Factored.out" ,
                         grm == "OctFullGRMGWAS" ~ "FullGrm",
                         grm == "OctLDGRMGWAS" ~ "LD",
                         grm == "OctLocoGRMGWAS" ~ "Loco",
                         grm == "OctNoGRMGWAS" ~ "None",
                         T ~ grm)
  ) %>% 
  mutate(trait = case_when(trait == "GIF" ~ "GIF",
                           trait == "LowestPvalue" ~ "LowestPvalue",
                           trait == "SNpcount" ~ "#ofHits")) %>% 
  filter(grm != "None",
         trait == "#ofHits") %>% 
  as.data.table(.)

fullprop2 = fullprop2 %>% 
  filter(chromosome %in% c("2L",  "3R"),
         inverted == "Inverted",
         grm %in% c( "Loco", "Factored.out")) 
%>%
  mutate(sig.pheno = ifelse(sig.pheno == T, "Inv traits", "Non-Inv traits"))
# saveRDS(fullprop2, "C:/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2025_objects/prop.beat.fig5")
g1colors = c("blue","red")

annotation_df <- data.frame(
  chromosome = c("2L","3R", "2L","3R"),
 # trait = c(rep("#ofHits",2),rep("GIF",2)),
  end = c(rep(2.1,2), rep(1.1,2)),
  start = c(rep(1.9,2),rep(0.9,2)),
  y = c(0.6, 0.45, 0.2, 0.22),
  label = c("**","NS.","NS.","NS. ")
)



j2 =   ggplot() +  
  ylab("Proportion Beating Permutations") +
  xlab("Methods") +
  # scale_color_manual(values = colors)+
  facet_grid(. ~ chromosome) +
  scale_color_manual( values = g1colors) +
  geom_errorbar(data = fullprop2, aes(
    x= grm,
    y=p,
    ymin=lowerci,
    ymax=upperci,
    color = sig.pheno),
    width = 0.1, position=position_dodge(width = 0.5), show.legend = F) +
  geom_point(data = fullprop2, aes(
    y= p,
    x=grm,
    color = sig.pheno),position=position_dodge(width = 0.5), show.legend = T) + theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  geom_vline(xintercept = 0.05, color = "grey", linetype = "dashed")+
  # geom_signif(
  #   data = annotation_df,
  #   aes(xmin = start, xmax = end, annotations = label, y_position = y),
  #   textsize = 5, vjust = -0.02,
  #   tip_length = 0.1,
  #   manual = TRUE
  # ) +
  coord_flip() + ylim(-0.1, 0.7)
j2


ggsave( ((g1 ) + plot_annotation(tag_levels = 'A') +
         plot_layout(guides = "collect") & theme(legend.position = 'bottom'))   , file = "C:/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2025_objects/images/gwasfig4.4.pdf", height = 6, width = 9) 

#####################
##Fig 4 statistics###
#####################
#use fishers exact test to compare proportions of signifcant pheontypes
#A compare between inv and std, and between loco and factored out. 
#B compare between sigtraits/nah, and between loco and factored out. 
#gifdata$phenotype = "GIF"
snpdata$phenotype = "SNP"
loppdata$phenotype = "lowP"
talldata = rbind( snpdata)

fixtraits = keeptraits %>% 
  rename("pheno"="phenotype")
talldata = talldata %>% 
  mutate(grm = case_when(grm == "OctcovarGRMGWAS" ~ "Factored.out" ,
                         grm == "OctFullGRMGWAS" ~ "FullGrm",
                         grm == "OctLDGRMGWAS" ~ "LD",
                         grm == "OctLocoGRMGWAS" ~ "Loco",
                         grm == "OctNoGRMGWAS" ~ "None",
                         T ~ grm)) %>% 
  filter(chromosome != 4) %>% 
  #rename("phenotype"= "pheno") %>% 
  merge(fixtraits, ., by = "pheno" )
  
ref.table = expand.grid(unique(talldata$grm), unique(talldata$chromosome), unique(talldata$pheno))
#gifdata comparisoins
g.out = foreach(f = c(1:dim(ref.table)[1]), .errorhandling = "remove")%do% {
  f = 1
  ref.info = ref.table[f,]
  dt = talldata[chromosome == unlist(ref.info[1,2])][grm == unlist(ref.info[1,1])][pheno == unlist(ref.info[1,3])]
   en =  dt[,
      list(TT=sum(dt == "inverted" & sig  == "Sig"),
           TF=sum(dt == "inverted" & sig  != "Sig"),
           FT=sum(dt != "inverted" & sig  == "Sig"),
           FF=sum(dt != "inverted" & sig  != "Sig"))]

  en[,or:=(TT/TF)/(FT/FF)]
  en[,fet.p:=apply(en, 1, function(x) fisher.test(matrix(as.numeric(c(x[1], x[2], x[3], x[4])), byrow=T, nrow=2))$p.value)]
  en$chromosome = unlist(ref.info[1,2])
  en$grm = unlist(ref.info[1,1])
  en$pheno = unlist(ref.info[1,3])
  en
}
g.out = rbindlist(g.out)
#I lowkey think we can just use the list feature for this. 
#signicant findings
gsig = g.out[fet.p < 0.05]

#next we want to check signicance of difference between loco and factored out, for every trait, for every inversion, for every chromosome
talldata = talldata %>% 
  mutate(inverted = gsub("Inverted","inverted",inverted),
         inverted = gsub("Non","non",inverted)) %>% 
  as.data.table(.)
en = talldata[,
  list(TT=sum(grm == "Loco" & sig == "Sig"),
       TF=sum(grm == "Loco" & sig == "Non-Sig"),
       FT=sum(grm == "Factored.out" & sig == "Sig"),
       FF=sum(grm == "Factored.out" & sig == "Non-Sig")),
  list(chromosome, inverted, phenotype)]
en[,or:=(TT/TF)/(FT/FF)]
en[,fet.p:=apply(en, 1, function(x) fisher.test(matrix(as.numeric(c(x[4], x[5], x[6], x[7])), byrow=T, nrow=2))$p.value)]
#wow that was worlds easier
sig.en = en[fet.p < 0.05]
#now we want to find the difference in suprassing perms between inv-associate and nah.
#we do it both for factored out and not.
#first- we need to bind inv associaton as a factor

head(talldata)
chr.ref = sig.phenos %>% 
  select(phenotype, inversion) %>% 
  #filter(inversion %in% c("In(2L)t","In(3R)Mo")) %>% 
  mutate(chromosome = case_when(inversion == "In(2L)t" ~ "2L", 
                                inversion == "In(2R)NS" ~ "2R", 
                                T ~ "3R"))
p.out = foreach(f = unique(chr.ref$chromosome), .combine = "rbind") %do% {
  #f = chr.ref$chromosome[1]
  chr.phenos = chr.ref %>% 
    filter(chromosome == f) 
  chr.phenos = unlist(chr.phenos$phenotype)
  delta.filter = talldata %>% 
    filter(chromosome == f) %>% 
    mutate(sig.pheno = pheno %in% chr.phenos)
}

#now we do more fischers exact test, across chromosome, summary stat, method, inversion group
en = p.out[,
              list(TT=sum(sig.pheno == T & sig == "Sig"),
                   TF=sum(sig.pheno == T & sig == "Non-Sig"),
                   FT=sum(sig.pheno == F & sig == "Sig"),
                   FF=sum(sig.pheno == F & sig == "Non-Sig")),
              list(chromosome, inverted, phenotype,grm)][inverted == "inverted"][grm %in% c("Factored.out","Loco")][chromosome %in% c("3R", "2L")]
en[,or:=(TT/TF)/(FT/FF)]
en[,fet.p:=apply(en, 1, function(x) fisher.test(matrix(as.numeric(c(x[5], x[6], x[7], x[8])), byrow=T, nrow=2))$p.value)]

sig.en.sig = en[fet.p < 0.05]

