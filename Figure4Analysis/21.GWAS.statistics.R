###run gif and number of sig snp analysis for each gwas output, within and without inversion on different chromosomes


### libraries
#library(GMMAT)
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel(4)
#library(ggplot2)
#get path to each gwas output
pheno.dir1 <- "/scratch/bal7cg/fall2023objects/OctcovarGRMGWAS/"
pheno.dir2 <- "/scratch/bal7cg/fall2023objects/OctFullGRMGWAS/"
pheno.dir3 <- "/scratch/bal7cg/fall2023objects/OctLDGRMGWAS/"
pheno.dir4 <- "/scratch/bal7cg/fall2023objects/OctLocoGRMGWAS/"
pheno.dir5 <- "/scratch/bal7cg/fall2023objects/OctNoGRMGWAS/"
inversions = readRDS("/scratch/bal7cg/inv.dt")
fixed3r = c("3R",11750567, 29031297,"total3r")

inversions = inversions %>% add_row(chrom = fixed3r[1], start = 11750567, stop = 29031297, invName = "fixed3")
#remove old rows
inversions = inversions[-c(3:5),]

#get paths to each gwas file
pheno.files <- list.files(path= c(pheno.dir1,pheno.dir2,pheno.dir3,pheno.dir4,pheno.dir5) , all.files=F, full.names=T, recursive=T )
pheno.filespaste = paste0("i=", pheno.files)

### define significant p values to look at
p <- 0.00001

o.statistics <- foreach(i = c(2:1311),.combine = "rbind", .errorhandling = "remove")%dopar%{
  #i = pheno.files[10]
  print(i)
  #load in score file
  gmmat = readRDS(i)
   #remove missing data
  gc(gmmat)
  
  #find genomic inflation
  
  chromout <- foreach(chr = unique(c("2L","2R","3L","3R")),.combine = "rbind") %do% {
    
   # chr = "2L"
    #filter down to a certain chromsomome
    gmmat.chr = gmmat[CHR == chr]
    #filter to within/without inversion, within first. 
    inv.breakpoint5 = inversions[chrom == chr]$start
    inv.breakpoint3 = inversions[chrom == chr]$stop
    inverted = gmmat.chr[POS <= inv.breakpoint3][POS >= inv.breakpoint5]
    #find genomic inflation
    chisq <- qchisq(1-inverted$PVAL,1)
    inv.GE = median(chisq)/qchisq(0.5,1)
    
    
    #find # of snps with pval < 0.001
    inv.sig.number = dim(inverted[PVAL < p]) [1]
    #find # of snps total
   inv.totalsnpnumber = dim(inverted)[1]
   #find for non-inverted region
   
   standard = gmmat.chr[POS > inv.breakpoint3 | POS < inv.breakpoint5]
   #find genomic inflation
   chisq <- qchisq(1-standard$PVAL,1)
   std.GE = median(chisq)/qchisq(0.5,1)
   
   
   #find # of snps with pval < 0.001
   std.sig.number = dim(standard[PVAL < p]) [1]
   #find # of snps total
   std.totalsnpnumber = dim(standard)[1]
    #put together summary statistics
    values = data.table (
      chromosome = chr ,
      inverted = c("inverted", "non-inverted"),
      GIF = c(inv.GE,std.GE),
      Sig.number = c(inv.sig.number,std.sig.number),
      Sig.threshold = p,
      Totalsnpnumber = c(inv.totalsnpnumber,std.totalsnpnumber),
      pheno = unlist(strsplit(i,"/",))[7],
      grm = unlist(strsplit(i,"/",))[5]
    )
    values
  }
  chromout
  
}
### change OUTPUTS and file names according to the real name
save(o.statistics,file  = "scratch/bal7cg/gwas.gif")



