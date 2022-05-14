################################################################
## FGFP MR: BMI --> Microbiome
## -- Pull data together --
## by: David Hughes
## date: March 16th 2022
################################################################
library(tidyverse)

################################
## Read in paramater file
################################
source("parameters/pfile.sh")
source("functions/knn_impute.R")
source("functions/rntransform.R")

################################
##
## Read in Source Data
##
################################
cat(paste0("1. Reading in the data\n"))
#################################
## 1) linker file
#################################
linker = read.table(linker_file, header = TRUE, as.is = TRUE, sep = "\t")

#################################
### 2) clinical data
#################################
clinical_data = read.table(clinical_data_file, header = TRUE, sep = "\t", as.is = TRUE)

## remove erroneous hip and waist circumference
w = which(colnames(clinical_data) %in% c("waist_circumference.1", "hip_circumference.1" ) )
if(length(w) > 0){ clinical_data = clinical_data[, -w] }

# alt_clinical_data = read.table(alt_clinical_data_file, header = TRUE, sep = "\t", as.is = TRUE)
# rownames(alt_clinical_data) = gsub("VDP.","VDP.0", rownames(alt_clinical_data))

#################################
### 3) GWASed microbiome data
#################################
gwased_mt_data = read.table(mt_gwased_traits_file, header = TRUE, sep = "\t", as.is = TRUE)


#################################
### 4) All microbial trait data
#################################
mt_data = read.table(microbiome_source_file, header = TRUE, sep = "\t", as.is = TRUE)

batch_data = read.table(lab_batch_file, header = TRUE, sep = "\t", as.is = TRUE)

#################################
### 5) MT phylogenetics
#################################
phylo_data = read.table(mt_phyogenetics_file, header = TRUE, sep = "\t", as.is = TRUE)


#################################
### 6) GRS estimates
#################################
grs_data = read.table(GRS_file, header = TRUE, sep = "\t", as.is = TRUE)


#################################
### 7) Metabolon Data
#################################
metabolon_data = read.table(metabolon_file, header = TRUE, sep = "\t", as.is = TRUE)
metabolon_feature_data = readr::read_tsv(metabolon_feature_file, show_col_types = FALSE)
metabolon_sample_data = readr::read_tsv(metabolon_sample_file, show_col_types = FALSE)
##
w = which(!metabolon_sample_data$SAMPLE_NAME %in% linker$VERC_id) 
metabolon_sample_data = metabolon_sample_data[-w,]
metabolon_data = metabolon_data[-w,]
##
m = match( rownames(metabolon_data) , linker$VERC_id )
rownames(metabolon_data) = linker$fgfp_id[m]
metabolon_sample_data$fgfp_id = linker$fgfp_id[m]


#################################
### 8) Restrict to overlapping samples
#################################
cat(paste0("2. Restricting and ordering by sample IDs\n"))
ids = sort( rownames(mt_data) )

## clinical_data
m = match(ids, rownames(clinical_data) )
clinical_data = clinical_data[m,]

# m = match(ids, rownames(alt_clinical_data) )
# alt_clinical_data = alt_clinical_data[m,]

## gwased_mt_data
m = match(ids, rownames(gwased_mt_data) )
gwased_mt_data = gwased_mt_data[m,]

## mt_data
m = match(ids, rownames(mt_data) )
mt_data = mt_data[m,]

## BATCH Data
m = match(ids, batch_data$SAMPLE_NAME )
batch_data = batch_data[m,]
rownames(batch_data) = batch_data$SAMPLE_NAME 
batch_data = batch_data[,-1]

## grs_data
x = rownames(grs_data); x = sapply( x, function(n){ strsplit( n, split = "_" )[[1]][1] } )
m = match(ids, x )
grs_data = grs_data[m,]
rownames(grs_data) = ids

## mt_data
m = match(ids, rownames(metabolon_data) )
metabolon_data = metabolon_data[m,]
## how many samples do not have metabolon data ?
sum( is.na(m) )
rownames( metabolon_data ) = ids

#################################
## X?) remove samples younger than 18
#################################
# w = which(clinical_data$age < 18)
# # length(w) ## 6
# if(length(w)>0){
#   clinical_data = clinical_data[-w, ]
#   gwased_mt_data = gwased_mt_data[-w, ]
#   mt_data = mt_data[-w, ]
#   grs_data = grs_data[-w, ]
#   metabolon_data = metabolon_data[-w, ]
#   batch_data = batch_data[-w,]
# }

#################################
## 9) Define WHRadjBMI
#################################
cat(paste0("3. Defining WHR and WHRadjBMI\n"))
# alt_clinical_data$whr = alt_clinical_data$Buikoptrek/ alt_clinical_data$Heupomtrek
##
clinical_data$WHR = clinical_data$waist_circumference / clinical_data$hip_circumference

fit = lm(WHR ~ BMI, data = clinical_data)
res = fit$residuals
m = match(rownames(clinical_data), names(res))
clinical_data$WHRadjBMI = res[m]


#################################
## 10) remove samples with no MT data
#################################
cat(paste0("4. Removing samlpes with NO microbiome data\n"))
missing = apply( mt_data[, c(3:505)] , 1, function(x){ sum( is.na(x) ) })
w = which(missing>0)  ## n = 2
##
if(length(w)>0){
  clinical_data = clinical_data[-w, ]
  gwased_mt_data = gwased_mt_data[-w, ]
  mt_data = mt_data[-w, ]
  grs_data = grs_data[-w, ]
  metabolon_data = metabolon_data[-w, ]
  batch_data = batch_data[-w,]
}

#################################
## 11) Pull the data together
#################################
cat(paste0("5. Pull everything together\n"))
## PULL EVERYTHING TOGETEHR
study_data = list()
study_data$clinical_data = clinical_data
study_data$gwased_mt_data = gwased_mt_data
study_data$mt_data = mt_data
study_data$batch_data = batch_data
study_data$grs_data = grs_data
study_data$phylo_data = phylo_data
study_data$metabolon_data = metabolon_data
study_data$metabolon_feature_data = metabolon_feature_data

#################################
## 12) SAVE data to Rdata file
#################################
cat(paste0("6. write to file\n"))
n = paste0( project_data_dir , "study_data_v1.Rdata")
save(study_data, file = n)


#################################
##
## 13) IMPUTATION
##
#################################
cat(paste0("7. Imputation\n"))
## how much missingness is there?
vars = c("age","sex","Weight","Height","BMI","hip_circumference","waist_circumference", "WHR", "WHRadjBMI")
missing_count = apply( clinical_data[, vars ], 2, function(x){sum(is.na(x))} )
# age                 sex              Weight              Height                 BMI   hip_circumference 
#  0                  30                  55                  56                  60                 153 
# waist_circumference                 WHR            WHRadjBMI 
#       151                           164                 171 

##########################
### load the class package
library(class)

#######################
## retain data with least 
## missingness
#######################
cat(paste0("7a. Limiting Metabolon, Microbiome, and Clincial variables to those useful\n"))
#######################
## METABOLITE DATA
#######################
## Identify metabolites with the least missingness
metabo_mis = apply(metabolon_data, 2, function(x) {sum(is.na(x))} )
w = which(metabo_mis == min(metabo_mis))
working_metabolon = metabolon_data[,w]
dim(working_metabolon) ## 412

#######################
## MICROBIOME DATA
#######################
## Identify microbial traits with the least missingness
mt_mis = apply(mt_data, 2, function(x) {sum(is.na(x))} )
w = which(mt_mis == min(mt_mis))
working_mt = mt_data[,w]
dim(working_mt) ## 829

## remove AB and PA
w = grep("PA", colnames(working_mt))
working_mt = working_mt[,-w]
dim(working_mt) ## 615

w = grep("AB", colnames(working_mt))
working_mt = working_mt[,-w]
dim(working_mt) ## 505

## Identify microbioal traits with the least zero
zero_prop = apply(working_mt, 2, function(x) { sum(x == 0, na.rm = TRUE)/length( na.omit(x) )   } )
w = which(zero_prop < 0.1 )
working_mt = working_mt[,w]
dim(working_mt) ## 413

## keep only those with a mean of 25
mean_est = apply(working_mt, 2, function(x) {  mean( as.numeric(x), na.rm = TRUE)   } )
w = which(mean_est > 40)
working_mt = working_mt[,w]
for(i in 1:ncol(working_mt)){ working_mt[,i] = as.numeric(working_mt[,i])}
## Data reduction of MT data
Cmat = cor(working_mt, method = "sp")
Dmat = as.dist(1 - abs(Cmat) )
nj = hclust(Dmat, method = "complete")
# plot(nj, hang = -1)
k = cutree(nj, h = 0.2)
ks = unique(k)
taxa_2_keep = sapply(ks, function(i){
  w = which(k == i)
  return(names(w)[1])
})
working_mt = working_mt[, taxa_2_keep]
## Evaluate tree cut
Cmat = cor(working_mt, method = "sp")
Dmat = as.dist(1 - abs(Cmat) )
nj = hclust(Dmat, method = "complete")
# plot(nj, hang = -1)

#######################
## CLINICAL DATA
#######################
clin_vars = c("Ratio_chol_tot_HDL.chol_mgmg","BP.sys", "BP.dia", "Vit_B12ngL", "Ureum_mgdL", "Trombocyten_1000mm3",
              "Triglyceriden_mgdL", "Totaal_Eiwit_g.L", "RBC_milj_mm3", "Neutrofielen_mm3", "Monocyten_mm3", "Monocyten_.",
              "MCV_fl", "MCHC_g_dL", "MCH_pg", "Lymfocyten_mm3", "Lymfocyten_.", "LDL.chol_gem_mgdL", "Insuline_nuchter_mUL",
              "Hemoglobine_gdL", "Hematocriet_.", "HDL.chol_mgdL", "HbA1c.", "HbA1c_mmolmol", "GPT_UL", "GOT_UL", "Gluc_nuchter_mg.dL", 
              "GammaGlobulines.", "GammaGlobulines_gL", "Gamma.GT_UL", "Ferritine_µg.L","Eosinofielen_mm3", "Eosinofielen_.", "e.GFR", 
              "CRP_mgL", "Creatinine_mgdL", "Cortisol_ochtend_mcgdL", "CK_UL", "Chol_totaal_mg.dL", "Bilirubine_tot_mgdL", "Bilirubine_indir_mgdL", 
              "Bilirubine_dir_mgdL", "Bez.na1_uur_mm", "BètaGlobulines.", "BètaGlobulines_gL", "Basofielen_mm3", "Basofielen_.", "Alfa2Globulines.",
              "Alfa2Globulines_gL", "Alfa1Globulines.", "Alfa1Globulines_gL", "Albumine.", "Albumine_gL")
mis = apply(clinical_data[, clin_vars], 2, function(x){ sum(is.na(x)) } )
w = names( which(mis<50) )
working_clin_data = clinical_data[,  w]

###########################
### IMPUTATION DATA SET
###########################
anthro_vars = c("age", "sex","Weight","Height","BMI","hip_circumference","waist_circumference", "WHR", "WHRadjBMI")
############
imputation_data = cbind(clinical_data[, anthro_vars],
              working_clin_data, 
              working_mt, 
              working_metabolon)
############
variables_2_test_for_imputation = colnames(imputation_data)[-c(1:9)]

############################
##
## *** Perfrom Imputation ***
##
###########################
vars = c("age","sex","Weight","Height","BMI","hip_circumference","waist_circumference", "WHR", "WHRadjBMI")

########################
### SEX IMPUTATION
########################
cat(paste0("7b. SEX Imputation\n"))
sex_imputation = knn_impute(wdata = imputation_data, 
           var_2_impute = "sex",
           r2_threshold = 0.1,
           variables_2_test_for_imputation = variables_2_test_for_imputation )
## sumstats
unlist( sex_imputation[2:4] )
## average imputation error
sex_imputation[5] ## 6.666%

## Add Imputation to data set
clinical_data$imputed_sex = clinical_data$sex
m = match( names(sex_imputation[[1]]) , rownames(clinical_data) )
clinical_data$imputed_sex[m] = as.character( sex_imputation[[1]] )

########################
### Weight IMPUTATION
########################
cat(paste0("7c. WEIGHT Imputation\n"))
imputation = knn_impute(wdata = imputation_data, 
                            var_2_impute = "Weight",
                            r2_threshold = 0.1,
                            variables_2_test_for_imputation = variables_2_test_for_imputation )
## sumstats
unlist( imputation[2:4] )
## average imputation error
imputation[5] ## 7.99 (95% CI 0.53-30.04)

## Add Imputation to data set
clinical_data$imputed_weight = clinical_data$Weight
m = match( names(imputation[[1]]) , rownames(clinical_data) )
clinical_data$imputed_weight[m] = as.numeric( as.character( imputation[[1]] ) )



########################
### Height IMPUTATION
########################
cat(paste0("7d. HEIGHT Imputation\n"))
imputation = knn_impute(wdata = imputation_data, 
                        var_2_impute = "Height",
                        r2_threshold = 0.1,
                        variables_2_test_for_imputation = variables_2_test_for_imputation )
## sumstats
unlist( imputation[2:4] )
## average imputation error
imputation[5] ## 5.15 (95% CI 0.16-16.74)

## Add Imputation to data set
clinical_data$imputed_height = clinical_data$Height
m = match( names(imputation[[1]]) , rownames(clinical_data) )
clinical_data$imputed_height[m] = as.numeric( as.character( imputation[[1]] ) )


########################
### BMI IMPUTATION
########################
cat(paste0("7e. BMI Imputation\n"))
imputation = knn_impute(wdata = imputation_data, 
                        var_2_impute = "BMI",
                        r2_threshold = 0.1,
                        variables_2_test_for_imputation = variables_2_test_for_imputation )
## sumstats
unlist( imputation[2:4] )
## average imputation error
imputation[5] ## 2.61 (95% CI 0.20-10.17)

## Add Imputation to data set
clinical_data$imputed_bmi = clinical_data$BMI
m = match( names(imputation[[1]]) , rownames(clinical_data) )
clinical_data$imputed_bmi[m] = as.numeric( as.character( imputation[[1]] ) )



########################
### hip_circumference IMPUTATION
########################
cat(paste0("7f. HIP Cir. Imputation\n"))
imputation = knn_impute(wdata = imputation_data, 
                        var_2_impute = "hip_circumference",
                        r2_threshold = 0.05,
                        variables_2_test_for_imputation = variables_2_test_for_imputation )
## sumstats
unlist( imputation[2:4] ) ## # of vars used = 24
## average imputation error
imputation[5] ## 5.82 (95% CI 0.16-23.50)

## Add Imputation to data set
clinical_data$imputed_hip_circumference = clinical_data$hip_circumference
m = match( names(imputation[[1]]) , rownames(clinical_data) )
clinical_data$imputed_hip_circumference[m] = as.numeric( as.character( imputation[[1]] ) )


########################
### waist_circumference IMPUTATION
########################
cat(paste0("7g. WAIST Cir. Imputation\n"))
imputation = knn_impute(wdata = imputation_data, 
                        var_2_impute = "waist_circumference",
                        r2_threshold = 0.1,
                        variables_2_test_for_imputation = variables_2_test_for_imputation )
## sumstats
unlist( imputation[2:4] ) ## # of vars used = 28
## average imputation error
imputation[5] ## 7.09 (95% CI 0.36-26.40)

## Add Imputation to data set
clinical_data$imputed_waist_circumference = clinical_data$waist_circumference
m = match( names(imputation[[1]]) , rownames(clinical_data) )
clinical_data$imputed_waist_circumference[m] = as.numeric( as.character( imputation[[1]] ) )


########################
### WHR IMPUTATION
########################
cat(paste0("7h. WHR Imputation\n"))
imputation = knn_impute(wdata = imputation_data, 
                        var_2_impute = "WHR",
                        r2_threshold = 0.1,
                        variables_2_test_for_imputation = variables_2_test_for_imputation )
## sumstats
unlist( imputation[2:4] ) ## # of vars used = 6
## average imputation error
imputation[5] ## 0.069 (95% CI 0.0028-0.26)

## Add Imputation to data set
clinical_data$imputed_WHR = clinical_data$WHR
m = match( names(imputation[[1]]) , rownames(clinical_data) )
clinical_data$imputed_WHR[m] = as.numeric( as.character( imputation[[1]] ) )


########################
### WHRadjBMI IMPUTATION
########################
cat(paste0("8. Est WHRadjBMI\n"))
fit = lm(imputed_WHR ~ imputed_bmi, data = clinical_data)
res = fit$residuals
m = match(rownames(clinical_data), names(res))
clinical_data$imputed_WHRadjBMI = res[m]





#################################
## 14) Pull the data together
#################################
cat(paste0("9. Pull Data Together\n"))

## PULL EVERYTHING TOGETEHR
study_data = list()
study_data$clinical_data = clinical_data
study_data$gwased_mt_data = gwased_mt_data
study_data$mt_data = mt_data
study_data$batch_data = batch_data
study_data$grs_data = grs_data
study_data$phylo_data = phylo_data
study_data$metabolon_data = metabolon_data
study_data$metabolon_feature_data = metabolon_feature_data

#################################
## 15) SAVE data to Rdata file
#################################
cat(paste0("10. Write to file\n"))
n = paste0( project_data_dir , "study_data_v1.1.Rdata")
save(study_data, file = n)


