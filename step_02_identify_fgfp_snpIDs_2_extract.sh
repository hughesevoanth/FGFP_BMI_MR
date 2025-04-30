################################
## FGFP MR
## by: David Hughes
## date: Sept 11 2024
################################

##################################
## Step 2
## -----------
## working in R
## read in the SNPSTATS files to 
## match rsids and or snpids (chr_pos and alleles)
##################################

## module add languages/R/4.4.1

library(tidyverse)
library(data.table)

## Source the paramater file with file paths
source("parameter_file.txt")

## Define the path to the snpstats files
snpstats_dir = paste0(DATADIR, "snpstats/")
## SNP stat files
files = list.files(snpstats_dir)

## READ in the Yengo and Pulit SNPs
pgs_snps = read.table( file = "/user/home/dh16508/fgfp/data/yengo_pulit_snps.txt", header= TRUE, sep = "\t", as.is = TRUE )
pgs_snps = pgs_snps |> arrange(CHR, POS)

## Now work through the SNP identification chr by chr
chrs = unique(pgs_snps$CHR)

new_pgs_snps = c()

###
for(chr in chrs){
	cat(paste0("\nnow processing ", chr, "\n"))

	if(chr < 10){
		f = paste0( snpstats_dir, "data_chr0", chr, ".snpstats" )
		working_snpstats = fread(f)
	} else {
		f = paste0( snpstats_dir, "data_chr", chr, ".snpstats" )
		working_snpstats = fread(f)
	}

	### limit the yengo pulit snps to those on the chromosome 'chr'
	temp_pgs_snps = pgs_snps |> filter(CHR == chr)
	### identify the subset of fgfp SNPs that match by position
	out1 = working_snpstats[position %in% unique(temp_pgs_snps$POS) ]
	## Reduce to columns to keep
	out1 = out1[, c(1:6,12,13,8,17 )]

	## match by position
	m = match(temp_pgs_snps$POS, out1$position)
	colnames(out1) = paste0("FGFP_", c("SNPID","RSID","CHR","POS","A1","A2","FREQ_A1","FREQ_A2","HWE_P","INFO") )
	temp_pgs_snps = cbind(temp_pgs_snps, out1[m,])

	## Do the RSIDs match ?
	temp_pgs_snps$RSID_match = ifelse(temp_pgs_snps$SNP == temp_pgs_snps$FGFP_RSID, 1, 0)

	## Do the ALLELEs match ?
	x = (temp_pgs_snps$Tested_Allele == temp_pgs_snps$FGFP_A1 & temp_pgs_snps$Other_Allele == temp_pgs_snps$FGFP_A2) | (temp_pgs_snps$Tested_Allele == temp_pgs_snps$FGFP_A2 & temp_pgs_snps$Other_Allele == temp_pgs_snps$FGFP_A1)
	temp_pgs_snps$ALLELE_match = ifelse(x, 1, 0)

	## Do the Allele Freq match ?
	x = (temp_pgs_snps$Tested_Allele == temp_pgs_snps$FGFP_A1 & abs(temp_pgs_snps$Freq_Tested_Allele_in_HRS - temp_pgs_snps$FGFP_FREQ_A1) < 0.15) | (temp_pgs_snps$Tested_Allele == temp_pgs_snps$FGFP_A2 & abs(temp_pgs_snps$Freq_Tested_Allele_in_HRS - temp_pgs_snps$FGFP_FREQ_A2) < 0.15)
	temp_pgs_snps$ALLELE_FREQ_match = ifelse(x, 1, 0)

	new_pgs_snps = rbind(new_pgs_snps, temp_pgs_snps)

}


write.table(new_pgs_snps, file = "yengo_pulit_snps_fgfp.txt", 
	row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

#################################
### Make a single file for all snps that will be used to generate any one GRS. 
### These will be used to run LD pruning. 
#################################
library(tidyverse)
mydata = read.table("yengo_pulit_snps_fgfp.txt", header = TRUE, sep = "\t", as.is = TRUE)
mydata = mydata |> arrange(PUB, TRAIT, CHR, POS )
traits = unique(mydata$TRAIT)

### YENGO BMI SNPs Primary and Secondary
temp = mydata |> filter(PUB == "yengo" & TRAIT == "bmi")
snps = temp |> dplyr::select(FGFP_RSID)
write.table(snps, file = "yengo_bmi_cojo_rsids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### YENGO BMI SNPs Primary ONLY
snps = temp |> filter(SECONDARY_SNP == 0) |> dplyr::select(FGFP_RSID)
write.table(snps, file = "yengo_bmi_primaryonly_rsids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### YENGO HEIGHT SNPs Primary and Secondary
temp = mydata |> filter(PUB == "yengo" & TRAIT == "height")
snps = temp |> dplyr::select(FGFP_RSID)
write.table(snps, file = "yengo_height_cojo_rsids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### YENGO HEIGHT SNPs Primary ONLY
snps = temp |> filter(SECONDARY_SNP == 0) |> dplyr::select(FGFP_RSID)
write.table(snps, file = "yengo_height_primaryonly_rsids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### BMI

### PULIT BMI SNPs INDEX SNPs
temp = mydata |> filter(PUB == "pulit" & TRAIT == "bmi")
snps = temp |> dplyr::select(FGFP_RSID)
write.table(snps, file = "pulit_bmi_index_rsids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### PULIT BMI MALES SNPs INDEX SNPs
temp = mydata |> filter(PUB == "pulit" & TRAIT == "bmi_males")
snps = temp |> dplyr::select(FGFP_RSID)
write.table(snps, file = "pulit_bmi_males_index_rsids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### PULIT BMI FEMALES SNPs INDEX SNPs
temp = mydata |> filter(PUB == "pulit" & TRAIT == "bmi_females")
snps = temp |> dplyr::select(FGFP_RSID)
dim(snps)
write.table(snps, file = "pulit_bmi_females_index_rsids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### WHR

### PULIT WHR SNPs INDEX SNPs
temp = mydata |> filter(PUB == "pulit" & TRAIT == "whr")
snps = temp |> dplyr::select(FGFP_RSID)
dim(snps) ## 316
write.table(snps, file = "pulit_whr_index_rsids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### PULIT WHR MALES SNPs INDEX SNPs
temp = mydata |> filter(PUB == "pulit" & TRAIT == "whr_males")
snps = temp |> dplyr::select(FGFP_RSID)
dim(snps) ## 79
write.table(snps, file = "pulit_whr_males_index_rsids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### PULIT WHR FEMALES SNPs INDEX SNPs
temp = mydata |> filter(PUB == "pulit" & TRAIT == "whr_females")
snps = temp |> dplyr::select(FGFP_RSID)
dim(snps) ## 203
write.table(snps, file = "pulit_whr_females_index_rsids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### whradjbmi

### PULIT whradjbmi SNPs INDEX SNPs
temp = mydata |> filter(PUB == "pulit" & TRAIT == "whradjbmi")
snps = temp |> dplyr::select(FGFP_RSID)
dim(snps) # 346
write.table(snps, file = "pulit_whradjbmi_index_rsids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### PULIT WHR MALES SNPs INDEX SNPs
temp = mydata |> filter(PUB == "pulit" & TRAIT == "whradjbmi_males")
snps = temp |> dplyr::select(FGFP_RSID)
dim(snps) ## 91
write.table(snps, file = "pulit_whradjbmi_males_index_rsids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

### PULIT whradjbmi FEMALES SNPs INDEX SNPs
temp = mydata |> filter(PUB == "pulit" & TRAIT == "whradjbmi_females")
snps = temp |> dplyr::select(FGFP_RSID)
dim(snps) ## 266
write.table(snps, file = "pulit_whradjbmi_females_index_rsids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)




