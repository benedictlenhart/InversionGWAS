###use linear modeling to find the association of cosmopoltian inversions with traits
#use a paired analysis to find mutations of the same frequency as inversions, and model them in the same way

# 1. make sure that the "same frequency mutations" also are being compared with the amount of dgrp lines as the corresponding inversion.
# if they have less lines , cut them.
# same amount, we're good.
# more lines, subsample down to the correct amount.
# 2. When determining which mutations are the same frequency, make sure that we are determining frequency of both -inversion and -mutation based on only the homozygotes- cut heterozygoes first. 
library(data.table)
library(foreach)
library(doParallel)

library(tidyverse)
library(SeqArray)
#cl = 4
registerDoParallel(10)
args = commandArgs(trailingOnly=TRUE)
#create job id from array
id = as.numeric(args[1])
id = 30
#load in inversion st and phenotype data
INPUT = "/project/berglandlab/Adam/fall2023objects/"#path of working directory

#load in reftable
reftable = readRDS(paste0(INPUT,"dgrpreftable"))

#load in our phenotype file 
phenotable = readRDS(paste0(INPUT, "sepphenolongform"))



#######################################
###make inv/same frequency models #####
#######################################



#merge in inversion info as factors
mergedata = merge(phenotable, reftable, by = "DGRP")
phenotype = unique(mergedata$fullpheno)[id]#487 total
pdata = mergedata %>% 
  
  filter(fullpheno == phenotype) %>% 
  select(-WolbachiaStatus_NA) %>% 
  pivot_longer(
    cols = c(Inversion_2L_t_NA,Inversion_2R_NS_NA ,Inversion_3R_K_NA, Inversion_3R_Mo_NA, Inversion_3R_P_NA),
    names_to = "Inversions",
    values_to = "Inv.st"
  ) %>% 
  na.omit(.) %>% 
  mutate(Inv.st = as.factor(Inv.st)) %>% 
  as.data.table(.)
#now we want to make a model for each inversion
inv.out = foreach(f = c("Inversion_2L_t_NA","Inversion_2R_NS_NA" ,"Inversion_3R_K_NA", "Inversion_3R_Mo_NA", "Inversion_3R_P_NA"))%dopar% {
  f = "Inversion_3R_P_NA"

  
  filtered.data = pdata[Inversions == f]
  #remove heterozygotes
  filtered.data = filtered.data[Inv.st != 1]
  samplenumber = dim(na.omit(filtered.data))[1]#exract out number of samples
  #find 
  
  #use reformulate to write formula
  form <- reformulate(termlabels = "Inv.st",response= "pheno.mean")
  testmodel = lm(form, data = filtered.data)
  
  #remove intercept and beta values
  inversiontype = unlist(testmodel$xlevels)
  beta = unlist(testmodel$coefficients)
  
  #create a null model and use anova to compare
  nullmodel= lm(pheno.mean ~ 1, filtered.data)
  anova = anova(testmodel, nullmodel)
  anova2 = anova(testmodel)
  #return a dataframe with phenotype, inversion, inversion type, beta, P value, and F statistic data
  o <- data.frame( 
    phenotype,
    inversion = f,
    P.value = unlist(anova[2,6]),
    F.stat = unlist(anova[2,5]),
    rsquared = summary(testmodel)$r.squared,
    adj.rsquared = summary(testmodel)$adj.r.squared
  )
  ###
  #now find
}
model.obs= rbindlist(inv.out)

#next step- find the frequency of each inversion within this specific phenotype's sample set, along with # of samples left after removing hets

#samplesize = length(unique(pdata$DGRP))
inv.count = pdata %>% 
  na.omit(.) %>% 
  filter(Inv.st != 1) %>% 
  group_by(Inversions) %>% 
  mutate(samplesize = n()) %>% 
  filter(Inv.st == 2) %>% 
  reframe(frequency = n()/ samplesize,
            samplesize = samplesize) %>% 
  distinct()
#this gives us the frequencies, and number of samples not including hets for each inversion
#here's the plan- iterated for loops (I know I know)
#first loop goes through each inversion, finds snps with similar frequency. second loop does modeling on each snp

inv2 = foreach(g = c(1:dim(inv.count)[1])) %do% {
  g = 5
  inv = unlist(inv.count[g,1])
  freq = unlist(inv.count[g,2])
  samplesize = unlist(inv.count[g,3])
  #we want to filter out hits that are in the inversion or near
  #in this version, we're making no snps within 2 mb of the inv breakpoits
  # I think alan said to do this including snps that could be within the inversion- lets just do that if we have time.
  #start filtering gds file 
  seqResetFilter(genofile)
  variant_info <- seqGetData(genofile, c("chromosome", "position"))
  chromosomes <- variant_info$chromosome
  positions <- variant_info$position
  if(inv == "Inversion_2L_t_NA") {
    chrom_filter <- chromosomes == "2L"
    position_filter <- positions < 225744 | positions > 15154180
    variants_to_keep <- which(chrom_filter & position_filter)
    seqSetFilter(genofile, variant.id = variants_to_keep)
  }else if ( inv == "Inversion_2R_NS_NA") {
    chrom_filter <- chromosomes == "2R"
    position_filter <- positions < 13391154 | positions > 22276334
    variants_to_keep <- which(chrom_filter & position_filter)
    seqSetFilter(genofile, variant.id = variants_to_keep)
    
  }else if ( inv == "Inversion_3R_K_NA") {
    chrom_filter <- chromosomes == "3R"
    position_filter <- positions < 9750567  | positions > 28140370
    variants_to_keep <- which(chrom_filter & position_filter)
    seqSetFilter(genofile, variant.id = variants_to_keep)
    
   
  }else if ( inv == "Inversion_3R_Mo_NA") {
    chrom_filter <- chromosomes == "3R"
    position_filter <- positions < 19406917  | positions > 31031297
    variants_to_keep <- which(chrom_filter & position_filter)
    seqSetFilter(genofile, variant.id = variants_to_keep)
    
  }else if ( inv == "Inversion_3R_RP_NA") {
    chrom_filter <- chromosomes == "3R"
    position_filter <- positions < 225744 | positions > 15154180
    variants_to_keep <- which(chrom_filter & position_filter)
    seqSetFilter(genofile, variant.id = variants_to_keep)
    
  }
  
 
  variant_sample_counts <- seqApply(genofile, "genotype", 
                                    function(geno) sum(is.na(geno)), 
                                    margin = "by.variant",
                                    as.is = "integer")#find number of non-missing samples for each variant
  #only include variants with atleast as many samples
  variants_to_keep <- which(variant_sample_counts >= samplesize)
  
  #now filter to only include variants with atleast as many samples as #samplecount
  seqSetFilter(genofile, variant.id = variants_to_keep)
  
  samples = seqGetData(genofile, "sample.id")
  dosage = seqGetData(genofile, "$dosage")
  #make a data object with samples and dosage
  dos = as.data.table(cbind(samples, dosage))
  dos = melt(dos, id.vars = c("samples"), variable.name = "variants", value.name = "dosage")
  #awesome, now we filter out all sample/mutant dosages that are missing or heterozygous
  dos = dos %>% 
    filter(dosage %in% c(0,2)) %>% 
    group_by(variants) %>% 
    mutate(sample.count = n()) %>% 
    filter(sample.count >= samplesize)#only include variants with atleast as many samples as the inversion
  #now filter to only include variants with similiar frequency
  keep = dos %>% 
    group_by(variants) %>% 
    filter(dosage == 0) %>% 
    mutate(mut.count = n()) %>% 
    reframe(mut.freq = mut.count / sample.count) %>% 
    filter(mut.freq >= (freq - 0.01),
           mut.freq <= (freq + 0.01),) %>% 
    unique()
  #merge in phenotype data
  phenotype.thinned = phenotable %>% 
    filter(fullpheno == phenotype)
  merge2 = dos %>% 
    merge(., keep, by= "variants") %>% 
    mutate(samples = gsub("line", "DGRP",samples)) %>% 
    rename("DGRP"="samples") %>% 
    merge(., phenotype.thinned, by= "DGRP")
  
  #once again filter down to only include variants with enough samples
  merge3 = merge2 %>% 
    na.omit() %>% 
    group_by(variants) %>% 
    mutate(new.count = n()) %>% 
    filter(new.count >= samplesize) %>% 
    as.data.table(.)
    
  #if we have greater than 500 variants left at this point, randomly pick 500. otherwise, use what we have
  if (length(unique(merge3$variants))[1]> 500) {
    print("toomany")
    variants2keep = sample(unique(merge3$variants), 500)
    merge4 = merge3 %>% 
      
      filter(variants %in% variants2keep)
  }else {
    variants2keep = unique(merge3$variants)
    merge4 = merge3
    print("toofew")
  }
  variants2keep = as.vector(variants2keep)
  #now create the second loop based on snp ids
  snp.out = foreach(s = c(1:length(variants2keep)), .errorhandling = "remove") %dopar% {
   #s = 1
    v = variants2keep[s]
    print(s)
    
    merge.dos = merge4 %>% 
      filter(variants == v)
    #now we need to downsample if the variant has excess samples
    if(length(unique(merge.dos$DGRP)) > samplesize) {
      samples2keep = sample(unique(merge.dos$DGRP), samplesize)
      merge.dos = merge.dos  %>% 
        filter(DGRP %in% samples2keep)
    }
    merge.dos$dosage = as.factor(merge.dos$dosage)
   
    testmodel = lm(pheno.mean ~ dosage, data = merge.dos)
    
    #remove intercept and beta values
    inversiontype = unlist(testmodel$xlevels)
    beta = unlist(testmodel$coefficients)
    
    #create a null model and use anova to compare
    nullmodel= lm(pheno.mean ~ 1, merge.dos)
    anova = anova(testmodel, nullmodel)
   
    #return a dataframe with phenotype, inversion, inversion type, beta, P value, and F statistic data
    o <- data.frame( 
      phenotype,
      inversion = inv,
      P.value = unlist(anova[2,6]),
      F.stat = unlist(anova[2,5]),
      rsquared = summary(testmodel)$r.squared,
      adj.rsquared = summary(testmodel)$adj.r.squared,
      snp.pos = s
    )
    
    o
  }
  freq.perms = rbindlist(snp.out)
  freq.perms
}
#bind from each inversion
freq.all = rbindlist(inv2)
#bind real and perm output into one file
model.obs$snp.pos = "inversion"
full.modeldata = rbind(freq.all, model.obs)


saveRDS(full.modeldata, paste0("/scratch/bal7cg/spring2025objects/freq.model/model", "-", phenotype))
