
library(data.table)
library(tidyverse)
library(foreach)
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Sep_2024_objects//")

#load in data on which phenotypes beat pod permutations
chr.ranked = readRDS("east.coast.chr.ranked")
chr.ranked.c = readRDS("chr.ranked")
#sig.phenos = readRDS("sigphenos")
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2025_objects//")
#load in data on which pheontypes we use, and which are significant. 
keeptraits = readRDS("keeptraits.fixed")
sig.phenos = readRDS("C:/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Feb_2025_objects/topmost.phenos")



chr.ranked$type = paste0("eastcoast.",chr.ranked$type)
chr.ranked.c$type = paste0("cville.", chr.ranked.c$type)

chr.ranked = rbind(chr.ranked, chr.ranked.c)

chr.keep = chr.ranked %>% 
  merge(keeptraits, ., by = "phenotype")

############################33
##compare to sig phenos
########################3
chr.ref = sig.phenos %>% 
  select(phenotype, inversion) %>% 
  mutate(chromosome = case_when(inversion == "In(2L)t" ~ "2L", 
                                inversion == "In(2R)NS" ~ "2R", 
                                T ~ "3R"))
#now create a foreach loop that filters delta and chr.ref by chr, then joins
chr.out = foreach(f = unique(chr.ref$chromosome), .combine = "rbind") %do% {
  #f = chr.ref$chromosome[1]
  chr.phenos = chr.ref %>% 
    filter(chromosome == f) 
  chr.phenos = unlist(chr.phenos$phenotype)
  delta.filter = chr.keep%>% 
    filter(chr == f) %>% 
    mutate(sig.pheno = phenotype %in% chr.phenos)
  delta.filter
}
#find proportion beating permutations
propsnp1 = chr.out %>% 
  #rename("phenotype"= "pheno") %>% 
  merge(keeptraits, ., by = "phenotype" ) %>% 
  group_by(sig.pheno, chr, inv.st, threshold, type, method) %>% 
  mutate(totalcount = n())
propsnp2 = propsnp1 %>% 
  #filter(sig == "Sig") %>% 
  #filter(avg.x > upperbound.y) %>% #we only want right tail
  group_by(sig.pheno, chr,threshold, type, method, totalcount,inv.st) %>% 
  summarize(pass.count = sum(sig == "Sig") )%>% 
  mutate(p = pass.count / totalcount,
         upperci = p + 1.96 * sqrt((p *(1-p)) / totalcount),
         lowerci = p - 1.96 * sqrt((p *(1-p)) / totalcount)) %>% 
  as.data.table(.)

propsnp2$type = gsub("enrichment", "", propsnp2$type)

propsnp2 = propsnp2 %>% 
  filter(chr %in% c("2L", "3R"),
         threshold == 500,
         inv.st == "Inverted",
         method %in% c( "Loco", "Factored.out"),
         type != "cville.annotation") %>% 
  mutate(sig.pheno = ifelse(sig.pheno == T, "Inv traits", "Non-Inv traits"),
         type = gsub("\\."," ", type),
         p = as.numeric(p)) 
  e2colors = c("blue","red")
  annotation_df <- data.frame(
    chr = c(rep("2L",8), rep("3R",8)),
    type = c(rep(
      c(rep("cville bayes ",2),rep("cville XtX ",2), rep("eastcoast bayes ",2), rep("eastcoast XtX ",2)),
      2)),
    end = c(rep(c(.9,1.9),8)),
    start = c(rep(c(1.1,2.1),8)), 
    y = c(rep(0.4, 4), rep(0.25, 6), 0.3, 0.5,rep(0.25, 4) ),
    label = c("NS.","*", rep(c("NS.", "NS. "), 4), "NS.","*", rep(c("NS.", "NS. "), 2))
  )
 


  
  e2 = ggplot() +  
    ylab("Proportion Exceeding Permutations") +
    xlab("Methods") +
    # scale_color_manual(values = colors)+
    facet_grid(chr ~ type) +
    scale_color_manual( values = e2colors) +
    geom_errorbar(data = propsnp2, aes(
      x= method,
      y=p,
      ymin=lowerci,
      ymax=upperci,
      color = sig.pheno),
      width = 0.1, position=position_dodge(width = 0.5), show.legend = F) +
    geom_point(data = propsnp2, aes(
      y= p,
      x=method,
      color = sig.pheno),position=position_dodge(width = 0.5), show.legend = T) + theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
    geom_vline(xintercept = 0.05, color = "grey", linetype = "dashed")+
    geom_signif(
      data = annotation_df,
      aes(xmin = start, xmax = end, annotations = label, y_position = y),
      textsize = 4, vjust = -0.02,
      tip_length = 0.15,
      manual = TRUE
    ) + ylim(-0.05, 0.55)
e2


ggsave(e2 , file = "C:/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2025_objects/images/gwasfig5.4.pdf", height = 6, width = 9) 

###############
##statistics##
##############



stats = chr.out %>% 
  #rename("phenotype"= "pheno") %>% 
  merge(keeptraits, ., by = "phenotype" ) %>%
  filter(threshold == 500) %>% 
  mutate(type = gsub("\\."," ", type),
         type = gsub("enrichment", "", type)) %>% 
  as.data.table(.)

#now we do more fischers exact test, across chromosome, summary stat, method, inversion group
en = stats[,
           list(TT=sum(sig.pheno == T & sig == "Sig"),
                TF=sum(sig.pheno == T & sig == "Non-Sig"),
                FT=sum(sig.pheno == F & sig == "Sig"),
                FF=sum(sig.pheno == F & sig == "Non-Sig")),
           list(chr, inv.st, type,method)] [inv.st == "Inverted"][method %in% c("Factored.out","Loco")][chr %in% c("3R", "2L")]
en[,or:=(TT/TF)/(FT/FF)]
en[,fet.p:=apply(en, 1, function(x) fisher.test(matrix(as.numeric(c(x[5], x[6], x[7], x[8])), byrow=T, nrow=2))$p.value)]

sig.en = en[fet.p < 0.05]
