library(lme4)
#library(aomisc)
library(data.table)
library(tidyverse)
library(scales)
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Sep_2023_objects/")
pdata = fread("all_phenotype_mean_data.tsv", header = T)

phenotypes = fread("phenotypes.tsv")

#load in metadata
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Sep_2024_objects/")
meta = fread("DGRP.meta.info.csv")
#simplify to one row per study
meta = meta %>% 
  select(!c("ph broad","narrow","V30")) %>% 
  distinct(.)
#we want columns for if inversions are mentioned, anovad, or factored out.
meta.grouped = meta %>% 
  rename("study_id"="Study ID") %>% 
  mutate(Factored_out = grepl( "nversion", Covariates),
         Anova = case_when(`In(2L)t ANOVA` == "Y" | `In(2R)NS ANOVA` == "Y" | `In(3R)P ANOVA` == "Y" | `In(3R)K ANOVA` == "Y"| `In(3R)Mo ANOVA` == "Y" ~ T,
                           T ~ F)) %>% 
  select(c(`Inversion?`,Factored_out, Anova, DOI, Study,study_id ))

table(meta.grouped$Factored_out)
#
#we want to refine phenotypes first
#start with only using fully integrated data
intdt = phenotypes[study_status == "integrated"]
#now we want to append the sex of each row to create sex-specific phenotypes
intdt = intdt %>% 
  select(c(phenotype_id, phenotype_name, nber_sex_female, nber_sex_male, nber_sex_na,study_id ))
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
#select down to just study_id, phenotype name, quantitity
phenotype.tall = meltdt %>% 
  rename("quantity_samples"="quantitiy") %>% 
  select(quantity_samples, study_id, phenotype_id,fullpheno)
#create merge data table
mergedt = merge(phenotype.tall, meta.grouped, by = "study_id")
#now merge in behavioral classications/filter to our final quality controlled traits
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Jan_2025_objects/")

keeptraits = readRDS("keeptraits.fixed")
merge2 = keeptraits %>%
  rename("fullpheno"="phenotype") %>% 
  mutate(meta.info = case_when(meta.info == "p"~ "Physiology",
                               meta.info == "sr"~ "Stress_Resistance",
                               meta.info == "lh"~ "Life_History",
                               meta.info == "m"~ "Morphology",
                               meta.info == "b"~ "Behavior")) %>%
  merge(mergedt, ., by = "fullpheno") %>% 
  as.data.table(.)

merge2 = merge2 %>% 
  relocate(Study, DOI, study_id, quantity_samples, `Inversion?`,Factored_out, Anova,phenotype_id, fullpheno, meta.info, V1)
#write new csv
write.csv(merge2, "Supp.meta.info.csv")

########################################
#how many pool studies are also in nunez 2024?
#############################################
#what are our doi's?
doi.pun = unique(merge2$DOI)
doi.pun = gsub("doi.org/", "", doi.pun)
#check which studies overlap with nunez 2024
nundt = fread("Table_S8_GENETICS-2023-306673.csv")
#get list of dois
doi.nun = unique(nundt$doi)
#remove doi script
doi.nun = gsub("doi:", "", doi.nun)
doi.nun = gsub("DOI:", "", doi.nun)
doi.nun = gsub("doi.", "", doi.nun)
doi.nun = gsub("DOI ", "", doi.nun)
doi.nun = gsub("org/", "", doi.nun)
doi.nun = gsub(" ", "", doi.nun)#remove spaces
#how many pool studies are also in nunez?

pool.in.n = doi.pun[(doi.pun %in% doi.nun)]
#how many nunez studies are also in pool?

nun.in.p = doi.nun[(doi.nun %in% doi.pun)]

###################################################
#check how many studies use factored out approach
######################################################
merge2 %>% 
  select(Study, Factored_out) %>% 
  distinct(.) %>% 
  filter(Factored_out == T)
