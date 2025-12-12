library(lme4)
#library(aomisc)
library(data.table)
library(tidyverse)
library(scales)
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Sep_2023_objects/")
#load in the data from the DGRPPool- we downloaded their phenotype data in fall of 2023
pdata = fread("all_phenotype_mean_data.tsv", header = T)

phenotypes = fread("phenotypes.tsv")



#we want to filter out phenotypes that are purely genomic features, or wolbachia, or have less then 75 lines.
#first lets remove traits we don't want. 

#we want to refine phenotypes firt
#start with only using fully integrated data
intdt = phenotypes[study_status == "integrated"]
#now we want to append the sex of each row to create sex-specific phenotypes
intdt = intdt %>% 
  select(c(phenotype_id, phenotype_name, nber_sex_female, nber_sex_male, nber_sex_na,study_id))
#use a melt to create rows for each sex
meltdt = melt(intdt, id.vars = c("phenotype_id", "phenotype_name", "study_id"), variable.name = "sex", value.name = "quantitiy" )
#for visual sake order by phenotype
meltdt = meltdt %>% 
  arrange(phenotype_id)
#filter out rows with no individuals, then create a new column with a pasted name
meltdt = meltdt %>% 
  filter(quantitiy > 0) %>% 
  mutate(phenotype.sex = case_when(sex == "nber_sex_female"~paste0(phenotype_id, "_F"),
                                   sex == "nber_sex_male"~ paste0(phenotype_id, "_M"),
                                   T ~ paste0(phenotype_id, "_NA") ) 
  ) %>% 
  mutate(fullpheno = case_when(sex == "nber_sex_female"~paste0(phenotype_name, "_F"),
                               sex == "nber_sex_male"~ paste0(phenotype_name, "_M"),
                               T ~ paste0(phenotype_name, "_NA") )     
  )


#lets check the distributions of samples use for each trait
ggplot(meltdt, aes(quantitiy)) + geom_histogram(bins = 15)
#very weird- most data is 125 or more, and then a big peak at 75. a definite tail sub 75, lets see how many
x = (meltdt[quantitiy < 75])
#67,, about 1/10 of the traits are sub 75. 

fullphenos = meltdt$fullpheno
wolb = fullphenos[grepl("olbachi", fullphenos) == T ]
cutices = fullphenos[grepl("Cuticul", fullphenos) == T]
variants = fullphenos[grepl("Variants", fullphenos) == T]
inversions = fullphenos[grepl("nversion", fullphenos) == T]
Vein = fullphenos[grepl("Vein", fullphenos) == T ]
vein = fullphenos[grepl("vein", fullphenos) == T ]
reads = fullphenos[grepl("reads", fullphenos) == T ]
Reads = fullphenos[grepl("Reads", fullphenos) == T ]
bp = fullphenos[grepl("Number_bp", fullphenos) == T ]
pelem = fullphenos[grepl("PElem", fullphenos) == T ]
genome = fullphenos[grepl("GenomeSize", fullphenos) == T ]
trans.elements = fullphenos[grepl("TEs", fullphenos) == T ]
binary = fullphenos[(fullphenos %in% c("Microbiota_OTU035_M", "Microbiota_OTU026_M")) ]
#bind names we want to remove
remove = c(cutices[-c(1:2)],wolb, variants, inversions, Vein, vein, "mn_Lifespan_F", reads,Reads,bp,genome,trans.elements, pelem, binary)#keeping 2 from cuticles, 155 total ro remove
#which phenotypes are left?
keep = subset(fullphenos, !(fullphenos %in% remove))
keepmelt = meltdt[fullphenos %in% keep]
#redo histogram
#filter to only include studies with atleast 75 samples
keepmelt = keepmelt[quantitiy >= 75]

#create long form of line averages
linemelt = melt(pdata, id.vars = c('DGRP'), variable.name = "phenotype.sex", value.name = "pheno.mean")
#filter to only include traits passing quality control
mergedata = merge(linemelt, keepmelt[,c(2,6)], by = "phenotype.sex")

#save long form phenotype data
saveRDS(mergedata, "sepphenolongform")
#note, this script reflects our final, most conservative quality control. Originally we included additional phenotypes-- scripts will use the "keeptraits" file to ensure the traits included in each final output are uniform with final quality control.
