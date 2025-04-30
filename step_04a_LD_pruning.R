################################################
## Compile the plink LD prun in out files
## and annotate the "yengo_pulit_snps_fgfp.txt" file
## indicating which SNPs are LD independent
## By: David A Hughes
## Date: Sep 13th 2024
################################################
library(tidyverse)
## Read in the SNP annotation file
mydata = read.table("../../data/gwas_snp_lists/20241111/yengo_pulit_snps_fgfp.txt", header = TRUE, sep = "\t", as.is = TRUE)

## Identify all of the ldprune files
datadir = "../../data/gwas_snp_lists/20241111/ldprune/"
files = list.files(datadir)

## Traits: 
traits = c("yengo_bmi_cojo.prune.in",
           "yengo_bmi_primaryonly.prune.in",
           "yengo_height_cojo.prune.in",
           "yengo_height_primaryonly.prune.in",
           "pulit_bmi.prune.in",
           "pulit_bmi_females.prune.in",
           "pulit_bmi_males.prune.in", 
           "pulit_whr.prune.in",
           "pulit_whr_females.prune.in",
           "pulit_whr_males.prune.in",
           "pulit_whradjbmi.prune.in",
           "pulit_whradjbmi_females.prune.in",
           "pulit_whradjbmi_males.prune.in")

snpcount = c()
## 
for(trait in traits){
  cojo = grep("cojo", trait)
  ## my trait name
  mytrait0 = gsub(".prune.in","",trait)
  mytrait = gsub("_cojo","",mytrait0)
  mytrait = gsub("_primaryonly","",mytrait)
  mypub = strsplit(mytrait, split = "_")[[1]][1]
  mytrait = paste( strsplit(mytrait, split = "_")[[1]][-1], collapse = "_" )
  
  ## Find the files
  fs = files[grep(trait, files)]
  ## read in each file and extract the rsids
  rsids = unlist( sapply(fs, function(x){
    read.table( paste0(datadir, x), header = FALSE)[,1] 
  }) )
  ## how many snps do I have in each ld pruned set?
  snpcount = c(snpcount, length(rsids))
  
  if(mypub == "yengo"){
    if(length(cojo) == 1){
      mydata = mydata |> mutate( !!mytrait0 := ifelse( TRAIT == mytrait & PUB == mypub, 1, 0) )  
    } else {
      mydata = mydata |> mutate( !!mytrait0 := ifelse( TRAIT == mytrait & PUB == mypub & SECONDARY_SNP == 0, 1, 0) )  
    }
  } else {
    mydata = mydata |> mutate( !!mytrait0 := ifelse( TRAIT == mytrait & PUB == mypub, 1, 0) )  
  }
  
  
  ## find the rows of data that match my trait AND rsid
  # w = which(mydata$TRAIT == mytrait & mydata$PUB == mypub & mydata$FGFP_RSID %in% rsids)
  mytrait0 = paste0("LD_ind_", mytrait0)
  mydata = mydata |> mutate( !!mytrait0 := ifelse( TRAIT == mytrait & PUB == mypub & FGFP_RSID %in% rsids, 1, 0) )
}

names(snpcount) = traits

## How many SNPs did I annotate for each trait?
w = grep("LD_ind_", colnames(mydata))
data.frame(source = snpcount,  anno = apply(mydata[,w], 2, sum) )
## Very good. perfect match !

## how many snps for all traits are there?
apply(mydata[,30:55], 2, sum) 


mydata |> filter(PUB == "yengo" & TRAIT == "bmi")

###################
## write to file
###################
write.table(mydata, file = "../../data/gwas_snp_lists/20241111/yengo_pulit_snps_fgfp_13092024.txt", 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

