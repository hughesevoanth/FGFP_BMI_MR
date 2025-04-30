################################
## FGFP MR
## by: David Hughes
## date: Sept 11 2024
################################

##################################
## Step 1
## -----------
## working in R
## extract rsids
##
##################################
library(tidyverse)

## Source the paramater file with file paths
source("parameters/snp_lists.sh")
files = ls()
## yengo file paths
w = grep("yengo", files )
yengo_files = files[w]
## pulit file paths
w = grep("pulit", files )
pulit_files = files[w]

#########################
## Merge Yengo Data
#########################
yb = read.table( eval(as.symbol(yengo_files[1])), header = TRUE, sep = "\t", as.is = TRUE)
yb = yb |> mutate(TRAIT = "bmi", .before = SNP)
yb = yb |> mutate(PUB = "yengo", .before = SNP)

yh = read.table( eval(as.symbol(yengo_files[2])), header = TRUE, sep = "\t", as.is = TRUE)
yh = yh |> mutate(TRAIT = "height", .before = SNP)
yh = yh |> mutate(PUB = "yengo", .before = SNP)

yengo_snps = rbind(yb, yh)
yengo_snps = yengo_snps |> mutate(SECONDARY_SNP = ifelse(P > 5e-8 , 1, 0 ) )

#########################
## Merge Pulit Data
#########################
## NOTE: allele "A1" in the Pulit data release is the TESTED allele. (Stated in supplementary data PDF of publication)

## BMI
## combined data
pb = read.table(eval(as.symbol("pulit_BMI")), header = TRUE, sep = " ", as.is = TRUE)
pb = pb[, c(1:9,11)]
colnames(pb) = colnames(yengo_snps)[3:12]
pb = pb |> mutate(TRAIT = "bmi", .before = SNP)
pb = pb |> mutate(PUB = "pulit", .before = SNP)

## males
pbm = read.table(eval(as.symbol("pulit_BMI_males")), header = TRUE, sep = " ", as.is = TRUE)
pbm = pbm[, c(1:5,20:23,25)]
colnames(pbm) = colnames(yengo_snps)[3:12]
pbm = pbm |> mutate(TRAIT = "bmi_males", .before = SNP)
pbm = pbm |> mutate(PUB = "pulit", .before = SNP)

## females
pbf = read.table(eval(as.symbol("pulit_BMI_females")), header = TRUE, sep = " ", as.is = TRUE)
pbf = pbf[, c(1:5,13:16,18)]
colnames(pbf) = colnames(yengo_snps)[3:12]
pbf = pbf |> mutate(TRAIT = "bmi_females", .before = SNP)
pbf = pbf |> mutate(PUB = "pulit", .before = SNP)

## WHR
## combined data
pw = read.table(eval(as.symbol("pulit_WHR")), header = TRUE, sep = " ", as.is = TRUE)
pw = pw[, c(1:9,11)]
colnames(pw) = colnames(yengo_snps)[3:12]
pw = pw |> mutate(TRAIT = "whr", .before = SNP)
pw = pw |> mutate(PUB = "pulit", .before = SNP)

## males
pwm = read.table(eval(as.symbol("pulit_WHR_males")), header = TRUE, sep = " ", as.is = TRUE)
pwm = pwm[, c(1:5,20:23,25)]
colnames(pwm) = colnames(yengo_snps)[3:12]
pwm = pwm |> mutate(TRAIT = "whr_males", .before = SNP)
pwm = pwm |> mutate(PUB = "pulit", .before = SNP)

## females
pwf = read.table(eval(as.symbol("pulit_WHR_females")), header = TRUE, sep = " ", as.is = TRUE)
pwf = pwf[, c(1:5,13:16,18)]
colnames(pwf) = colnames(yengo_snps)[3:12]
pwf = pwf |> mutate(TRAIT = "whr_females", .before = SNP)
pwf = pwf |> mutate(PUB = "pulit", .before = SNP)

## WHRadjBMI
## combined data
pwb = read.table(eval(as.symbol("pulit_whradjbmi")), header = TRUE, sep = " ", as.is = TRUE)
pwb = pwb[, c(1:9,11)]
colnames(pwb) = colnames(yengo_snps)[3:12]
pwb = pwb |> mutate(TRAIT = "whradjbmi", .before = SNP)
pwb = pwb |> mutate(PUB = "pulit", .before = SNP)

## males
pwbm = read.table(eval(as.symbol("pulit_whradjbmi_males")), header = TRUE, sep = " ", as.is = TRUE)
pwbm = pwbm[, c(1:5,20:23,25)]
colnames(pwbm) = colnames(yengo_snps)[3:12]
pwbm = pwbm |> mutate(TRAIT = "whradjbmi_males", .before = SNP)
pwbm = pwbm |> mutate(PUB = "pulit", .before = SNP)

## females
pwbf = read.table(eval(as.symbol("pulit_whradjbmi_females")), header = TRUE, sep = " ", as.is = TRUE)
pwbf = pwbf[, c(1:5,13:16,18)]
colnames(pwbf) = colnames(yengo_snps)[3:12]
pwbf = pwbf |> mutate(TRAIT = "whradjbmi_females", .before = SNP)
pwbf = pwbf |> mutate(PUB = "pulit", .before = SNP)


#### Pull all pulit files together
pulit_snps = rbind(pb, pbm, pbf, pw, pwm, pwf, pwb, pwbm, pwbf)
pulit_snps$BETA_COJO = NA
pulit_snps$SE_COJO = NA
pulit_snps$SE_COJO = NA
pulit_snps$P_COJO = NA
pulit_snps$SECONDARY_SNP = NA

## edit the pulit rsids snps
extracted_rsid = sapply(pulit_snps$SNP, function(x){strsplit(x , split = ":" )[[1]][1] })
pulit_snps$SNP = extracted_rsid


### All Yengo and Pulit SNPs together
anthro_snps = rbind(yengo_snps, pulit_snps)
write.table(anthro_snps, file = "yengo_pulit_snps.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE )



