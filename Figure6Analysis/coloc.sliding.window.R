#do a sliding window coloc analysis on pca-gwas outputs
#goal: load in the two gwas. then iterature through both of them using windows
#at each window, use color to find single causal varients. 
library(data.table)
#library(tidyverse)
library(foreach)
library(doParallel)
library(coloc)
library(matrixStats)

registerDoParallel(10)

#load data
#goal- regather gwas and report number of sig snps and gif accross all chromsoomes. 
pheno.dir2 <- "/project/berglandlab/Adam/fall2023objects/factored_out.pca.gwas2"

pheno.dir5 <- "/project/berglandlab/Adam/fall2023objects/loco.pca.gwas2"

# inversions = readRDS("/scratch/bal7cg/inv.dt")
# 
# pheno.files <- list.files(path=c(pheno.dir2, pheno.dir5), all.files=T, full.names=T, recursive=T)
#choose files we want
#we want 10 sets of files, dim1/2 for each inversion, for each method. 
chrom = c("In(2L)t","In(2R)NS" ,"In(3R)Mo","In(3R)K")
ref.table = data.frame(
  inversion = rep(chrom,1),
  method = rep(c(pheno.dir2, pheno.dir5), each = 4)
)


#create a frame foreach that takes them through each phenotype
toppheno = foreach(f = c(1:8), .combine = "rbind") %do% {
   #f = 1
  #use reftable to get the files we need
  refslice = ref.table[f,]
  inv = unlist(refslice[1,1])
 
  path = unlist(refslice[1,2])
  if(path == "/project/berglandlab/Adam/fall2023objects/factored_out.pca.gwas2"){
  pheno1 = readRDS(paste0(path, "/", inv, "-Dim.1.txt .RDS"))
  pheno2 = readRDS(paste0(path, "/", inv, "-Dim.2.txt .RDS"))
  }else {
    pheno1.2L = readRDS(paste0(path, "/", inv, "-Dim.1#2L"))
    pheno1.2R = readRDS(paste0(path, "/", inv, "-Dim.1#2R"))
    pheno1.3L = readRDS(paste0(path, "/", inv, "-Dim.1#3L"))
    pheno1.3R = readRDS(paste0(path, "/", inv, "-Dim.1#3R"))
    pheno1 = rbind(pheno1.2L, pheno1.2R, pheno1.3L, pheno1.3R)
    
    pheno2.2L = readRDS(paste0(path, "/", inv, "-Dim.2#2L"))
    pheno2.2R = readRDS(paste0(path, "/", inv, "-Dim.2#2R"))
    pheno2.3L = readRDS(paste0(path, "/", inv, "-Dim.2#3L"))
    pheno2.3R = readRDS(paste0(path, "/", inv, "-Dim.2#3R"))
    pheno2 = rbind(pheno2.2L, pheno2.2R, pheno2.3L, pheno2.3R)
  }
 
  
  # win.bp <- 1e6
  # step.bp <- 5e5
  # #create a table of windows using the gwas output file
  # tmp <- pheno1[CHR == chromosome]
  # wins = data.table(chr=chromosome,
  #            start=seq(from=min(tmp$POS), to=max(tmp$POS)-win.bp, by=step.bp),
  #            end=seq(from=min(tmp$POS), to=max(tmp$POS)-win.bp, by=step.bp) + win.bp)
  # wins[,i:=1:dim(wins)[1]]

    
  win.bp <- 1e4
  step.bp <- 5e3

  wins <- foreach(chr.i=c("2L","2R","3L","3R"), .combine="rbind", .errorhandling="remove")%do%{
    #chr.i = "2L"
    tmp <- pheno1[CHR == chr.i]
    data.table(chr=chr.i,
               start=seq(from=min(tmp$POS), to=max(tmp$POS)-win.bp, by=step.bp),
               end=seq(from=min(tmp$POS), to=max(tmp$POS)-win.bp, by=step.bp) + win.bp)
  }
  wins[,i:=1:dim(wins)[1]]
  dim(wins)  
  cwin<- foreach(win.i=1:dim(wins)[1], .combine = "rbind", .errorhandling = "remove")%dopar%{
    #if(win.i%%100==0) message(win.i)
    print(win.i)
   #  win.i <- 9
    #w
    # tmp <- gwas1[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]
    #not sure what Alan's doing so we're just going to keep it simple
    #use our window data table to slice out a specific region of the genome within the gwas
    g1slice = pheno1[CHR == wins[win.i]$chr][POS >= wins[win.i]$start][POS <= wins[win.i]$end]
    
    g2slice = pheno2[CHR == wins[win.i]$chr][POS >= wins[win.i]$start][POS <= wins[win.i]$end]
    #next, use this gwas to create our coloc- appropriate data set
    #need a list with components (betas, varbeta, snp(names), position, type, maf
    g1list = list(
      beta = g1slice$SCORE,
      varbeta = g1slice$VAR,
      snp = g1slice$SNP,
      postition = g1slice$POS,
      type = "quant",
      MAF = g1slice$AF,
      N = dim(g1slice)[1]
    )
    g2list = list(
      beta = g2slice$SCORE,
      varbeta = g2slice$VAR,
      snp = g2slice$SNP,
      postition = g2slice$POS,
      type = "quant",
      MAF = g2slice$AF,
      N = dim(g2slice)[1]
    )
    
    #check out what a potential causal locus might be just for g1
    #my.res <- finemap.abf(dataset=g1list)
    #hmm- many likely snps. linkage maybe?
    my.res <- coloc.abf(dataset1=g1list,
                        dataset2=g2list)
    #we want to return the H4 probabilty, and the position of interesting snps
    h4 = my.res$summary[6]
    h1 = my.res$summary[3]
    h2 = my.res$summary[4]
    snp.names = (subset(my.res$results,SNP.PP.H4>0.01))$snp
    snp.props = (subset(my.res$results,SNP.PP.H4>0.01))$SNP.PP.H4
    max.prop = max(snp.props)
    #for now, just do top snp
    if(is.null(snp.names) == T) {
      o = data.table(
        chr = wins[win.i]$chr,
        win.i,
        h4,
        h1,
        h2,
        snp.names = NA,
        snp.props = NA
      )
    } else {
      o = data.table(
        chr = wins[win.i]$chr,
        win.i,
        h4,
        h1,
        h2,
        snp.names = (subset(my.res$results,SNP.PP.H4>0.01))$snp,
        snp.props = (subset(my.res$results,SNP.PP.H4>0.01))$SNP.PP.H4
      )
    }
  
    o[,start:=wins[win.i]$start]
    o[,end:=wins[win.i]$end]
    o
  }
  cwin$inv = inv

  cwin$method = path
  cwin
}

saveRDS(toppheno,"/scratch/bal7cg/fall2024objects/colocpca.1.3.fullgenome")
