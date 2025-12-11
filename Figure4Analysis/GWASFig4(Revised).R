#Alan's figure creation- new gwas figure

### libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(foreach)
library(doMC)

fullout = readRDS("/project/berglandlab/Adam/etoh.summary.RDS")
gwas = readRDS("/project/berglandlab/Adam/etoh.g2")
fullprop2 = readRDS("/project/berglandlab/Adam/fall2023objects/gwasfig4.data")

colors = c("grey","cyan")
#create boxplot
g1 = ggplot()+
  geom_boxplot(data = fullout[perm.st == "Permuted"], aes(x = grm, y = pcount), color = "grey") +
  geom_jitter(data = fullout[perm.st == "Permuted"], aes(x = grm, y = pcount), color = "grey", alpha = 0.3)+
  geom_point(data = fullout[perm.st == "Observed"], aes(x = grm, y = pcount), color = "cyan", size = 10, shape = "diamond")+
  theme_bw(base_size = 15)+
  ylab("Number of Highly Significant SNPs") + xlab("GRM Method")+
  ylim(0,140)


g2colors = c("orange1","navyblue")
g2 = ggplot(data=gwas[PVAL<.05][pheno%in%c("mn_AlcoholTolerance_M")], aes(x=POS, y=-log10(PVAL), color=grm)) + 
  geom_vline(xintercept=14617051) + 
  geom_vline(xintercept=5943643) + 
  geom_vline(xintercept=5950033) +
  geom_vline(xintercept=2166548, linetype="dashed") +
  geom_vline(xintercept=13139098, linetype="dashed") +
  geom_point(alpha=.5) +
  scale_color_manual( values = g2colors) +
  facet_grid(pheno~.) + 
  scale_x_continuous(n.breaks = 3) + 
  theme_bw(base_size = 15)
#bring in original figure 4

library(ggsignif)
fullprop2 = readRDS("/project/berglandlab/Adam/fall2023objects/gwasfig4.data")
g3colors = c("blue","red")


annotation_df <- data.frame(
  chromosome = c(rep(c("2L","2R", "3L", "3R"), 3)),
  #trait = c(rep("#ofHits",4)),
  end = c(rep(2,4),rep(3,4),rep(4,4)),
  start = c(rep(1,12)),
  y = c(rep(0.15, 4),rep(0.3, 4),rep(0.45, 4)),
  label = c(rep("NS./NS. ", 4), rep("NS./NS.", 4), rep("***/***", 4))
)
g3 = ggplot() +  
  ylab("Proportion exceeding permutations") +
  xlab("Methods") +
  # scale_color_manual(values = colors)+
  facet_grid(. ~ chromosome) +
  scale_color_manual( values = g3colors) +
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
g3
library(patchwork)
ggsave( ((g3 / (g1 + g2) ) + plot_annotation(tag_levels = 'A') +
           plot_layout(guides = "collect") & theme(legend.position = 'bottom'))   ,
        file = "gwasfig4.5.png", height = 9, width = 9)



