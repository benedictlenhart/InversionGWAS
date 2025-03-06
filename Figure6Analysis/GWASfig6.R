#make a set of graphs that shows the h4 value comparing dim1/dim2 amongst each window. could also consider showing the number of potentially co-significant snps
library(data.table)
#library(tidyverse)
library(foreach)
library(doParallel)
library(tidyverse)
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Feb_2025_objects/")
cd = readRDS("colocpca.1.3.fullgenome")
cd = readRDS("colocpca.1.4.fullgenome")
cd = rbindlist(cd[-c(5,10)])#3rp failed due to so few data
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Nov_2024_objects/")
inv.dt = readRDS("inv.dt")
inv.dt = inv.dt %>% 
  mutate(invName = case_when(invName == "2Lt"~"In(2L)t",
                             invName == "2RNS"~"In(2R)NS",
                             invName == "3RK"~"In(3R)K",
                             invName == "3RMo"~"In(3R)Mo",
                             invName == "3RP"~"In(3R)P"))

colnames(inv.dt)[1] = "chr"
colnames(inv.dt)[4] = "inv"
# inv.dt1 = inv.dt[inv %in% c("In(2L)t","In(3R)Mo")]
# #inv.dt = inv.dt[chr != "3L"]
# #make mimic datasets for other 3r inversions
# inv2 = inv.dt[inv != "In(3R)Mo"]
# inv2$inv = gsub("In\\(3R\\)K","In(3R)Mo", inv2$inv)
# inv2 = inv2[inv == "In(3R)Mo"]
# 
# inv3 = inv.dt[inv != "In(3R)Mo"]
# inv3$inv = gsub("In\\(3R\\)P","In(3R)Mo", inv3$inv)
# inv3 = inv3[inv == "In(3R)Mo"]
inv.2L = inv.dt[1,]
inv.3R = inv.dt[chr == "3R"]
inv.3R$height = c(0.4, 0.8, 0.6)
groupdata = cd %>% 
  mutate(method = ifelse(grepl("factored",method) == T, "Factored_out", "LOCO")) %>% 
 # filter (inv %in% c("2L_t","3R_Mo")) %>% 
  group_by(chr,inv, win.i, h4,h1,h2, start, end,method) %>% 
  summarise(number.snps = n()) %>% 
  #we want to melt the h's
  pivot_longer(cols = c(h1, h2, h4), names_to = "Principal component", values_to = "Enrichment" ) %>% 
  mutate(`Principal component` = case_when(`Principal component` == "h1" ~ "PC1",
                                           `Principal component` == "h2" ~ "PC2",
                                           T ~ "PC1 & PC2")) %>% 
  #mutate(inv = gsub("_","",inv)) %>% 
  as.data.table(.)
l.dt = data.frame(
  label = c("In(2L)t"),
  x = c(0.145e7),
  y = c(.65),
  inv = c("In(2L)t"),
  chr = c("2L")
)
r.dt = data.frame(
  label = c("In(3R)K","In(3R)P","In(3R)Mo"),
  x = c( 1.4e7, 1.9e7,2.78e7),
  y = c( 0.52, 0.7, 0.6),
  inv = c( "In(3R)Mo","In(3R)Mo","In(3R)Mo"),
  chr = c("3R","3R","3R")
)

g1= ggplot() +
  geom_rect(data=inv.2L, aes(xmin=start, xmax=stop, ymin=0, ymax=.8), color="black", fill = "grey", alpha=.2) +
  # geom_rect(data=inv.3R, aes(xmin=start, xmax=stop, ymin=0, ymax= height), color="black", fill = "grey", alpha=.2) +
  geom_point(data= groupdata[method == "Factored_out"][chr == "2L"][inv == "In(2L)t"],
             aes(x=round((start + end)/2), y= Enrichment), color="red", size = 1)  +
  facet_grid(`Principal component`~., scales = "free") +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text( face = "bold",angle = 45,hjust = 1, vjust = 1, size = 12))+
  geom_text(data = l.dt,aes( x = x, y = y, label = label), size = 3) +
  scale_x_continuous(n.breaks = 4) + 
  xlab("Position")+ ylab("Propability of causality") +ylim(0,0.81)
g1
g2= ggplot() +
  # geom_rect(data=inv.2L, aes(xmin=start, xmax=stop, ymin=0, ymax=.8), color="black", fill = "red", alpha=.2) +
  geom_rect(data=inv.3R, aes(xmin=start, xmax=stop, ymin=0, ymax= height), color="black", fill = "grey", alpha=.2) +
  geom_point(data= groupdata[method == "Factored_out"][chr == "3R"][inv == "In(3R)Mo"],
             aes(x=round((start + end)/2), y= Enrichment), color="red", size = 1)  +
  facet_grid(`Principal component`~., scales = "free") +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text( face = "bold",angle = 45,hjust = 1, vjust = 1, size = 12))+
  geom_text(data = r.dt,aes( x = x, y = y, label = label), size = 3) +
  scale_x_continuous(n.breaks = 4) + 
  xlab("Position")+ ylab("Propability of causality") +ylim(0,0.81)
g1
g2
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2025_objects/images/")
ggsave((g1 + g2 + 
       plot_annotation(
         tag_levels = c("A"))) ,
       filename = "Fig6fac.pdf", width = 10, height = 7 )
#where are these peaks?
groupdata[method == "LOCO"][inv %in% c("2Lt","In(3R)Mo")][h4 > 0.4]
