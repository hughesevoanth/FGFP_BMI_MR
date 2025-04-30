###################################################
## FGFP BMI-Microbiome MR
## by: David Hughes
## date: Sep 13 2024
###################################################

######################################
##
## (I) load data
##
######################################
library(tidyverse)

## Read in the SNP annotation file
mydata = read.table("../../data/gwas_snp_lists/20241111/yengo_pulit_snps_fgfp_13092024.txt", header = TRUE, sep = "\t", as.is = TRUE)

## Read in my function
source("functions/grs_prepare_dosage_data.R")
source("functions/build_grs.R")

###################
## (Ia) READ IN Dosage data
###################
dosage_file = "../../data/geno/yengo_pulit_snps.dosage"
dosage_data = read.table( dosage_file , header = TRUE, as.is = TRUE)
dosage_data = as.data.frame(dosage_data)


######################################
##
## (II) Estimate the GRS for Yengo BMI
##  
######################################
## Primary and Secondary SNPs using the cojo estimates
temp = mydata |> filter(yengo_bmi_cojo == 1 & FGFP_INFO >= 0.3 & FGFP_HWE_P > 1e-8)

yengo_bmi_cojo = build_grs( grs_data = temp,
                  dos_data = dosage_data,
                  grs_rsids_colname = "FGFP_RSID" ,
                  grs_EA_colname = "Tested_Allele" ,
                  grs_AA_colname = "Other_Allele" ,
                  grs_EAF_colname = "Freq_Tested_Allele_in_HRS" ,
                  grs_beta_colname = "BETA_COJO"  )

## Primary SNPs 
temp = mydata |> filter(yengo_bmi_primaryonly == 1 & FGFP_INFO >= 0.3 & FGFP_HWE_P > 1e-8)

## Primary SNPs using the initial gwas estimates
yengo_bmi_primary = build_grs( grs_data = temp,
                                 dos_data = dosage_data,
                                 grs_rsids_colname = "FGFP_RSID" ,
                                 grs_EA_colname = "Tested_Allele" ,
                                 grs_AA_colname = "Other_Allele" ,
                                 grs_EAF_colname = "Freq_Tested_Allele_in_HRS" ,
                                 grs_beta_colname = "BETA"  )

plot(yengo_bmi_primary, yengo_bmi_cojo)
cor.test(yengo_bmi_primary, yengo_bmi_cojo)


######################################
##
## (II) Estimate the GRS for Yengo Height
##  
######################################
## Primary and Secondary SNPs using the cojo estimates
temp = mydata |> filter(yengo_height_cojo == 1 & FGFP_INFO >= 0.3 & FGFP_HWE_P > 1e-8)

yengo_height_cojo = build_grs( grs_data = temp,
                            dos_data = dosage_data,
                            grs_rsids_colname = "FGFP_RSID" ,
                            grs_EA_colname = "Tested_Allele" ,
                            grs_AA_colname = "Other_Allele" ,
                            grs_EAF_colname = "Freq_Tested_Allele_in_HRS" ,
                            grs_beta_colname = "BETA_COJO"  )

## Primary SNPs 
temp = mydata |> filter(yengo_height_primaryonly == 1 & FGFP_INFO >= 0.3 & FGFP_HWE_P > 1e-8)

## Primary SNPs using the initial gwas estimates
yengo_height_primary = build_grs( grs_data = temp,
                               dos_data = dosage_data,
                               grs_rsids_colname = "FGFP_RSID" ,
                               grs_EA_colname = "Tested_Allele" ,
                               grs_AA_colname = "Other_Allele" ,
                               grs_EAF_colname = "Freq_Tested_Allele_in_HRS" ,
                               grs_beta_colname = "BETA"  )

plot(yengo_height_primary, yengo_height_cojo)
cor.test(yengo_height_primary, yengo_height_cojo)



######################################
##
## (II) Estimate the GRS for Pulit BMI
##  
######################################
## Index SNPs
temp = mydata |> filter(pulit_bmi == 1 & FGFP_INFO >= 0.3 & FGFP_HWE_P > 1e-8)

pulit_bmi = build_grs( grs_data = temp,
                               dos_data = dosage_data,
                               grs_rsids_colname = "FGFP_RSID" ,
                               grs_EA_colname = "Tested_Allele" ,
                               grs_AA_colname = "Other_Allele" ,
                               grs_EAF_colname = "Freq_Tested_Allele_in_HRS" ,
                               grs_beta_colname = "BETA"  )

## Index SNPs 
temp = mydata |> filter(pulit_bmi_females == 1 & FGFP_INFO >= 0.3 & FGFP_HWE_P > 1e-8)

pulit_bmi_females = build_grs( grs_data = temp,
                                  dos_data = dosage_data,
                                  grs_rsids_colname = "FGFP_RSID" ,
                                  grs_EA_colname = "Tested_Allele" ,
                                  grs_AA_colname = "Other_Allele" ,
                                  grs_EAF_colname = "Freq_Tested_Allele_in_HRS" ,
                                  grs_beta_colname = "BETA"  )

## Index SNPs
temp = mydata |> filter(pulit_bmi_males == 1 & FGFP_INFO >= 0.3 & FGFP_HWE_P > 1e-8)

pulit_bmi_males = build_grs( grs_data = temp,
                               dos_data = dosage_data,
                               grs_rsids_colname = "FGFP_RSID" ,
                               grs_EA_colname = "Tested_Allele" ,
                               grs_AA_colname = "Other_Allele" ,
                               grs_EAF_colname = "Freq_Tested_Allele_in_HRS" ,
                               grs_beta_colname = "BETA"  )


plot(pulit_bmi, pulit_bmi_females)
cor.test(pulit_bmi, pulit_bmi_females)
cor.test(pulit_bmi, pulit_bmi_males)
##
cor.test(pulit_bmi, yengo_bmi_primary)
cor.test(pulit_bmi, yengo_bmi_cojo)



######################################
##
## (II) Estimate the GRS for Pulit WHR
##  
######################################
## Index SNPs
temp = mydata |> filter(pulit_whr == 1 & FGFP_INFO >= 0.3 & FGFP_HWE_P > 1e-8)

pulit_whr = build_grs( grs_data = temp,
                       dos_data = dosage_data,
                       grs_rsids_colname = "FGFP_RSID" ,
                       grs_EA_colname = "Tested_Allele" ,
                       grs_AA_colname = "Other_Allele" ,
                       grs_EAF_colname = "Freq_Tested_Allele_in_HRS" ,
                       grs_beta_colname = "BETA"  )

## Index SNPs 
temp = mydata |> filter(pulit_whr_females == 1 & FGFP_INFO >= 0.3 & FGFP_HWE_P > 1e-8)

pulit_whr_females = build_grs( grs_data = temp,
                               dos_data = dosage_data,
                               grs_rsids_colname = "FGFP_RSID" ,
                               grs_EA_colname = "Tested_Allele" ,
                               grs_AA_colname = "Other_Allele" ,
                               grs_EAF_colname = "Freq_Tested_Allele_in_HRS" ,
                               grs_beta_colname = "BETA"  )

## Index SNPs
temp = mydata |> filter(pulit_whr_males == 1 & FGFP_INFO >= 0.3 & FGFP_HWE_P > 1e-8)

pulit_whr_males = build_grs( grs_data = temp,
                             dos_data = dosage_data,
                             grs_rsids_colname = "FGFP_RSID" ,
                             grs_EA_colname = "Tested_Allele" ,
                             grs_AA_colname = "Other_Allele" ,
                             grs_EAF_colname = "Freq_Tested_Allele_in_HRS" ,
                             grs_beta_colname = "BETA"  )
###
cor.test(pulit_whr, pulit_whr_males)
cor.test(pulit_whr, pulit_whr_females)

######################################
##
## (II) Estimate the GRS for Pulit pulit_whradjbmi
##  
######################################
## Index SNPs
temp = mydata |> filter(pulit_whradjbmi == 1 & FGFP_INFO >= 0.3 & FGFP_HWE_P > 1e-8)

pulit_whradjbmi = build_grs( grs_data = temp,
                       dos_data = dosage_data,
                       grs_rsids_colname = "FGFP_RSID" ,
                       grs_EA_colname = "Tested_Allele" ,
                       grs_AA_colname = "Other_Allele" ,
                       grs_EAF_colname = "Freq_Tested_Allele_in_HRS" ,
                       grs_beta_colname = "BETA"  )

## Index SNPs 
temp = mydata |> filter(pulit_whradjbmi_females == 1 & FGFP_INFO >= 0.3 & FGFP_HWE_P > 1e-8)

pulit_whradjbmi_females = build_grs( grs_data = temp,
                               dos_data = dosage_data,
                               grs_rsids_colname = "FGFP_RSID" ,
                               grs_EA_colname = "Tested_Allele" ,
                               grs_AA_colname = "Other_Allele" ,
                               grs_EAF_colname = "Freq_Tested_Allele_in_HRS" ,
                               grs_beta_colname = "BETA"  )

## Index SNPs
temp = mydata |> filter(pulit_whradjbmi_males == 1 & FGFP_INFO >= 0.3 & FGFP_HWE_P > 1e-8)

pulit_whradjbmi_males = build_grs( grs_data = temp,
                             dos_data = dosage_data,
                             grs_rsids_colname = "FGFP_RSID" ,
                             grs_EA_colname = "Tested_Allele" ,
                             grs_AA_colname = "Other_Allele" ,
                             grs_EAF_colname = "Freq_Tested_Allele_in_HRS" ,
                             grs_beta_colname = "BETA"  )
###
cor.test(pulit_whradjbmi, pulit_whradjbmi_males)
cor.test(pulit_whradjbmi, pulit_whradjbmi_females)


grs_ests = data.frame(
  ids = names(yengo_height_primary),
  yengo_bmi_cojo = yengo_bmi_cojo,
  yengo_bmi_primary = yengo_bmi_primary, 
  yengo_height_cojo = yengo_height_cojo,
  yengo_height_primary = yengo_height_primary,
  pulit_bmi = pulit_bmi, 
  pulit_bmi_females = pulit_bmi_females, 
  pulit_bmi_males = pulit_bmi_males,
  pulit_whr = pulit_whr,
  pulit_whr_females = pulit_whr_females,
  pulit_whr_males = pulit_whr_males,
  pulit_whradjbmi = pulit_whradjbmi,
  pulit_whradjbmi_females = pulit_whradjbmi_females,
  pulit_whradjbmi_males = pulit_whradjbmi_males
)

write.table(grs_ests, file = "../../data/GRS/estimated_GRSs_20241113.txt", 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)



