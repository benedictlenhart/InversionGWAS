# create pca and loadings plots
library(data.table)
library(tidyverse)
library(foreach)
library(ggsignif)

wd =  paste0("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Sep_2024_objects/" )
scree.var = readRDS("C:/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Sep_2024_objects/scree.var")
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2025_objects/")
loadings = readRDS( "pca_loading.data2")
pca = readRDS( "pca_dim.data2")


loadings.small = loadings %>% 
  filter(inv %in% c("In(2L)t", "In(3R)Mo")) %>% 
  select(-c(V2, mu3, se3)) %>% 
  pivot_longer(cols = c(mu1, mu2), names_to = "PC", values_to = "mu") %>% 
  pivot_longer(cols = c(se1, se2), names_to = "PCse", values_to = "se") %>% 
  mutate(PC = gsub("mu","PC",PC)) %>% 
  as.data.table(.)
#there's got to be a better way to do this, but for now clear mismatches
loadings.small1 = loadings.small[PC == "PC1" & PCse == "se1"] 
loadings.small2 = loadings.small[PC == "PC2" & PCse == "se2"] 
loadings.small = rbind(loadings.small1, loadings.small2)
#make scree plot
p0colors = c("brown4","pink1")
p0 = scree.var %>% 
  filter(Inversion %in% c("In(2L)t","In(3R)Mo")) %>% 
  ggplot() +
  geom_line(aes(x = PC, y = eigen, group = Inversion, color = Inversion))+
  geom_point(aes(x = PC, y = eigen, group = Inversion, color = Inversion)) +
  theme_bw(base_size = 15) +
  scale_color_manual(values = p0colors) +
  ylab("Variance explained") + xlab("Principal Component")

p1colors = c("blue","red")
annotation_df <- data.frame(
  inv = c("In(2L)t","In(2L)t", "In(3R)Mo", "In(3R)Mo"),
  end = c(1.1, 2.1,1.1, 2.1),
  start = c(.9, 1.9,.9, 1.9),
  y = c(2, 2, 2, 4),
  label = c("***", "NS.","***", "*")
)
p1<- ggplot() + 
  geom_errorbar(data = loadings.small, aes(color=as.factor(gt.name), ymin=mu-1.96*se, ymax=mu+1.96*se, x = PC), width=.1, position = position_dodge(width= 0.3)) +
  geom_point(data=loadings.small, aes(color=as.factor(gt.name), y=mu, x = PC), position = position_dodge(width= 0.3)) +
  facet_grid(.~inv)+
  scale_color_manual(values = p1colors) +
  theme_bw(base_size = 15)  + ylim(-6,5) +
  xlab(NULL) +
  ylab("")+
  theme(legend.title =element_blank()) +
  geom_signif(
    data = annotation_df,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    textsize = 7, vjust = -0.2,
    tip_length = 0.01,
    manual = TRUE
  )


p1


lim <- 1.5
small = pca%>% 
  filter(inv %in% c( "In(3R)Mo")) %>% 
  mutate(pheno = gsub("completeObs.", "", pheno))
pR = ggplot(data=small) +
  coord_equal() +
  geom_segment(aes(x=0, y=0, xend=Dim.1, yend=Dim.2), arrow=arrow(length=unit(0.1,"cm")), size = 0.3) +
 
  xlim(-lim,lim) + ylim(-lim,lim) +
  xlab("PC1 (25.8 of Var)")+
  ylab("PC2 (14.8 of Var)")+
  theme_bw() +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))

small = pca%>% 
  filter(inv %in% c( "In(2L)t")) %>% 
  mutate(pheno = gsub("completeObs.", "", pheno))
pL = ggplot(data=small) +
  coord_equal() +
  geom_segment(aes(x=0, y=0, xend=Dim.1, yend=Dim.2), arrow=arrow(length=unit(0.1,"cm")), size = 0.3) +

  xlim(-lim,lim) + ylim(-lim,lim) +
  xlab("PC1 (21.2 of Var)")+
  ylab("PC2 (12.1 of Var)")+
  theme_bw() +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
pL
pR
ggsave(p3, file = "pcae2L_t.pdf", width = 3.61, height = 3.37)
# Saving 3.61 x 3.37 in image for 3rmo
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2025_objects/images/")
# ggsave( ((p3.3 | p3.3) /  p1 / p2) | g3  +plot_annotation(tag_levels = 'A'), file = "gwasfigpca.3.pdf", height = 14, width = 10) 
ggsave( ((p0 |p1) /(pL | pR) + plot_layout(guides = "collect") & theme(legend.position = 'bottom'))   +plot_annotation(tag_levels = 'A'), file = "gwasfig2.3.pdf", height = 12, width = 12) 

#peel through the traits, helps with identification of which traits map onto each quadrant of pc space
top.3r = small %>% 
  filter(y >0) %>% 
  arrange(x)
bottom.3r = small %>% 
  filter(y <0) %>% 
  arrange(x)

top.2L = small %>% 
  filter(Dim.2 > 0,
         Dim.1 < 0) %>% 
  select(c(Dim.1, Dim.2, pheno)) %>% 
  arrange(Dim.1)
bottom.2L = small %>% 
  filter(Dim.2 < 0) %>% 
  arrange(x)
inv = small %>% 
  filter(Dim.1 > -0.25,
         Dim.1 < 0.1,
         Dim.2 > 0)

std.3r = small %>% 
  filter(Dim.1 > 0.5)
