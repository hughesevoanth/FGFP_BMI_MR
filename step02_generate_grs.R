###################################################
## FGFP BMI-Microbiome MR
## by: David Hughes
## date: March 14 2022
###################################################

######################################
##
## (I) load data
##
######################################
library(tidyverse)

source("functions/grs_prepare_dosage_data.R")

###################
## (Ia) READ IN GWAS SNP info
###################
source("parameters/snp_lists.sh")
files = ls()
## yengo file paths
w = grep("yengo", files )
yengo_files = files[w]
## pruit file paths
w = grep("pruit", files )
pruit_files = files[w]

yengo_data = lapply(yengo_files, function(file){
	d = read.table(eval(as.symbol(file)), header = TRUE, sep = "\t", as.is = TRUE)
	d = as.data.frame(d)
	## filter for PRIMARY SNPs
	d = d %>% filter(P<=1e-8 & P_COJO<=1e-8)
	return(d)
})
names(yengo_data) = yengo_files

pruit_data = lapply(pruit_files, function(file){
	d = read.table(eval(as.symbol(file)), header = TRUE, sep = " ", as.is = TRUE)
	## remove alleles from SNP names
	d$SNP = sapply(d$SNP, function(x){ strsplit(x, split = ":")[[1]][1] })
	d = as.data.frame(d)
	return(d)
})
names(pruit_data) = pruit_files

###################
## (Ib) READ IN Dosage data
###################
dosage_data = read.table( dosage_file , header = TRUE, as.is = TRUE)
dosage_data = as.data.frame(dosage_data)
snpstats = read.table( snpstats_file , header = TRUE, as.is = TRUE, skip = 10)
snpstats = as.data.frame(snpstats)
# sample_data = read.table( sample_file , header = TRUE, as.is = TRUE)


######################################
##
## (II) Estimate the GRS for Yengo
##	
######################################
GRS_ests_Yengo = list()
for(n in names(yengo_data) ){
	GRS_ests_Yengo[[ n ]] = grs_prepare_dosage_data( 
		dosage_data = dosage_data, 
		snpstats_data = snpstats, 
		grs_data  = yengo_data[[ n ]] , 
		grs_data_SNP_col_name = "SNP" , 
		grs_data_EffectAllele_col_name = "Tested_Allele",
		grs_data_AltAllele_col_name = "Other_Allele",
		grs_data_EffectAllele_FREQ_col_name = "Freq_Tested_Allele_in_HRS",
		grs_data_BETA_col_name = "BETA",
		name_of_grs_being_processed = n )
}



######################################
##
## (III) Estimate the GRS for Pruit
##	
######################################
GRS_ests_Pruit = list()
for(n in names(pruit_data)[c(1,4,7)] ){
	GRS_ests_Pruit[[ n ]] = grs_prepare_dosage_data( 
		dosage_data = dosage_data, 
		snpstats_data = snpstats, 
		grs_data  = pruit_data[[ n ]] , 
		grs_data_SNP_col_name = "SNP" , 
		grs_data_EffectAllele_col_name = "A1.combined",
		grs_data_AltAllele_col_name = "A2.combined",
		grs_data_EffectAllele_FREQ_col_name = "frqA1.combined",
		grs_data_BETA_col_name = "beta.combined",
		name_of_grs_being_processed = n )
}

## FEMALES
GRS_FEMALES_ests_Pruit = list()
for(n in names(pruit_data)[c(2,5,8)] ){
	GRS_FEMALES_ests_Pruit[[ n ]] = grs_prepare_dosage_data( 
		dosage_data = dosage_data, 
		snpstats_data = snpstats, 
		grs_data  = pruit_data[[ n ]] , 
		grs_data_SNP_col_name = "SNP" , 
		grs_data_EffectAllele_col_name = "A1.combined",
		grs_data_AltAllele_col_name = "A2.combined",
		grs_data_EffectAllele_FREQ_col_name = "frqA1.combined",
		grs_data_BETA_col_name = "beta.females",
		name_of_grs_being_processed = n )
}

GRS_MALES_ests_Pruit = list()
for(n in names(pruit_data)[c(3,6,9)] ){
	GRS_MALES_ests_Pruit[[ n ]] = grs_prepare_dosage_data( 
		dosage_data = dosage_data, 
		snpstats_data = snpstats, 
		grs_data  = pruit_data[[ n ]] , 
		grs_data_SNP_col_name = "SNP" , 
		grs_data_EffectAllele_col_name = "A1.combined",
		grs_data_AltAllele_col_name = "A2.combined",
		grs_data_EffectAllele_FREQ_col_name = "frqA1.combined",
		grs_data_BETA_col_name = "beta.males",
		name_of_grs_being_processed = n )
}

#########################
## GRS summary statistics
##
#########################
a = sapply(GRS_ests_Yengo, function(x){x[[2]]})
b = sapply(GRS_ests_Pruit, function(x){x[[2]]})
d = sapply(GRS_FEMALES_ests_Pruit, function(x){x[[2]]})
e = sapply(GRS_MALES_ests_Pruit, function(x){x[[2]]})

grs_sumstats = cbind(a,b,d,e)
rownames(grs_sumstats) = rownames(GRS_ests_Yengo[[1]][[2]] )

write.table(grs_sumstats, file = "GRS_SumStats.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)


#########################
## GRS Estimates
##
#########################
a = c(); for(i in 1:length(GRS_ests_Yengo) ){ a = cbind(a,GRS_ests_Yengo[[i]][[1]]) }
b = c(); for(i in 1:length(GRS_ests_Pruit) ){ b = cbind(b,GRS_ests_Pruit[[i]][[1]]) }
d = c(); for(i in 1:length(GRS_FEMALES_ests_Pruit) ){ d = cbind(d,GRS_FEMALES_ests_Pruit[[i]][[1]]) }
e = c(); for(i in 1:length(GRS_MALES_ests_Pruit) ){ e = cbind(e,GRS_MALES_ests_Pruit[[i]][[1]]) }

grs_est = cbind(a,b,d,e)

write.table(grs_est, file = "GRS_Estimates.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

