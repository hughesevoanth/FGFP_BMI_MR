################################################################
## FGFP MR: BMI --> Microbiome
## -- Pull data together --
## by: David Hughes
## date: Sept 17th 2024
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
clinical_data = readxl::read_xlsx(clinical_data_file, sheet = 1)
nrow(clinical_data) ## 2998

## edit the sample IDs
clinical_data$vdp_ids = gsub("VDP.","VDP.0", clinical_data$vdp_ids )

## PCs and clinical blood work
JW_clinical_data = read.table(clinical_data_file_vJW, header = TRUE, sep = "\t", as.is = TRUE)
nrow(JW_clinical_data) ## 3133

pcs = JW_clinical_data[, paste0("PC", 1:10)]

blood_traits = JW_clinical_data[, c("Albumine_gL","Albumine.","Alfa1Globulines_gL",
                                    "Alfa1Globulines.","Alfa2Globulines_gL",
                                    "Basofielen_.","Basofielen_mm3",
                                    "BètaGlobulines_gL","BètaGlobulines.",
                                    "Bez.na1_uur_mm","Bilirubine_dir_mgdL","Bilirubine_indir_mgdL",
                                    "Bilirubine_tot_mgdL","Chol_totaal_mg.dL","CK_UL",
                                    "Cortisol_ochtend_mcgdL", "Creatinine_mgdL", "CRP_mgL",
                                    "e.GFR", "Eosinofielen_.", "Eosinofielen_mm3",
                                    "Ferritine_µg.L","Gamma.GT_UL","GammaGlobulines_gL",
                                    "GammaGlobulines.","Gluc_nuchter_mg.dL","GOT_UL",
                                    "GPT_UL","HbA1c_mmolmol","HbA1c.",
                                    "HDL.chol_mgdL","Hematocriet_.","Hemoglobine_gdL",
                                    "Insuline_nuchter_mUL", "LDL.chol_gem_mgdL","Lymfocyten_.",
                                    "Lymfocyten_mm3","MCH_pg", "MCHC_g_dL",
                                    "Monocyten_.","Monocyten_mm3",
                                    "Neutrofiele_segmenten_.","Neutrofielen_mm3","Ratio_chol_tot_HDL.chol_mgmg",
                                    "RBC_milj_mm3","Totaal_Eiwit_g.L", "Triglyceriden_mgdL",
                                    "Trombocyten_1000mm3","Ureum_mgdL","Urinezuur_mgdL",
                                    "Vit_B12ngL", "WBC_mm3") ]

# Cmat = cor(blood_traits, method = "spearman", use = "pairwise.complete.obs")
# Dmat = as.dist( 1 - abs(Cmat) )
# # Dmat = as.dist(1 - Cmat)
# nj = hclust(Dmat, method = "complete")
# plot(nj, hang = -1)
# k = cutree(nj, h = 0.5)


#################################
### 2a) drug data
#################################
## Level 2 ATC codes
drug_data = readxl::read_xlsx(drugs_data_level2_file, sheet = 1)
nrow(drug_data) ## 3049

## edit the sample IDs
drug_data$vdp_clean_id = gsub("VDP.","VDP.0", drug_data$vdp_clean_id )

## Level 3 ATC codes
drug_data_l3 = readxl::read_xlsx(drugs_data_level3_file, sheet = 1)
nrow(drug_data_l3) ## 3049

## edit the sample IDs
drug_data_l3$vdp_clean_id = gsub("VDP.","VDP.0", drug_data_l3$vdp_clean_id )


#################################
### 3) GWASed microbiome data
#################################
gwased_mt_data = read.table(mt_gwased_traits_file, header = TRUE, sep = "\t", as.is = TRUE)
nrow(gwased_mt_data) ## 2259

#################################
### 4) All microbial trait data
#################################
mt_data = read.table(microbiome_source_file, header = TRUE, sep = "\t", as.is = TRUE)
nrow(mt_data) ## 2259

batch_data = read.table(lab_batch_file, header = TRUE, sep = "\t", as.is = TRUE)
nrow(batch_data) ## 2482

#################################
### 5) MT phylogenetics
#################################
phylo_data = read.table(mt_phyogenetics_file, header = TRUE, sep = "\t", as.is = TRUE)
nrow(phylo_data)

#################################
### 6) GRS estimates
#################################
grs_data = read.table(GRS_file, header = TRUE, sep = "\t", as.is = TRUE)
nrow(grs_data) ## 2242

new_grs_data = read.table(new_GRS_file, header = TRUE, sep = "\t", as.is = TRUE)
nrow(new_grs_data) ## 2548

#################################
### 7) Metabolon Data
#################################
metabolon_data = read.table(metabolon_file, header = TRUE, sep = "\t", as.is = TRUE)
nrow(metabolon_data) ## 2993
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
ids = sort( rownames(mt_data)[rownames(mt_data) %in% clinical_data$vdp_ids] )
ids = ids[ids %in% rownames(metabolon_data) ]
length(ids) ## 2239
#### These individuals are not in the clincial data release: "VDP.02036" "VDP.02205" "VDP.03452"

## clinical_data
m = match(ids, clinical_data$vdp_ids )
clinical_data = clinical_data[m,]

m = match(ids, rownames(pcs))
pcs = pcs[m,]

m = match(ids, rownames(blood_traits))
blood_traits = blood_traits[m,]

## drug data
m = match(ids, drug_data$vdp_clean_id)
drug_data = drug_data[m,]

m = match(ids, drug_data_l3$vdp_clean_id)
drug_data_l3 = drug_data_l3[m,]

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

## new_grs_data
x = new_grs_data$ids; x = sapply( x, function(n){ strsplit( n, split = "_" )[[1]][1] } )
m = match(ids, x )
new_grs_data = new_grs_data[m,]
rownames(new_grs_data) = ids

## The new GRS estimates correlate negatively.
## Let's flip the distribution of values.
new_grs_data[,-1] = apply(new_grs_data[,-1], 2, function(x){ 
  max(x, na.rm = TRUE) - (x - min(x,  na.rm = TRUE))
  })
cor.test(new_grs_data$yengo_bmi_primary, clinical_data$BMI_online)$estimate ## 0.19
cor.test(new_grs_data$yengo_bmi_primary, grs_data$yengo_BMI_wGRS_info_hwe_filtered)$estimate ## 0.9893013

## metabolon_data
m = match(ids, rownames(metabolon_data) )
metabolon_data = metabolon_data[m,]
## how many samples do not have metabolon data ?
sum( is.na(m) )
rownames( metabolon_data ) = ids

#################################
## X?) remove samples younger than 18
#################################
w = which(clinical_data$age_mq < 18)
length(w) ## 9
if(length(w)>0){
  clinical_data = clinical_data[-w, ]
  pcs = pcs[-w,]
  blood_traits = blood_traits[-w, ]
  drug_data = drug_data[-w, ]
  drug_data_l3 = drug_data_l3[-w, ]
  gwased_mt_data = gwased_mt_data[-w, ]
  mt_data = mt_data[-w, ]
  grs_data = grs_data[-w, ]
  new_grs_data = new_grs_data[-w, ]
  metabolon_data = metabolon_data[-w, ]
  batch_data = batch_data[-w,]
}


#################################
## 9) Limit Drug data to just the antibiotic codes
##    J01, A07A, S01A, S02A
##      and to diarrhea (A07) and constipation (A06)
#################################
## A07 = 96
## A06 = 33

## table(drug_data$J01)
## J01 = 624

## There are no S02A users
## A07A = 2
## S01A = 24
## S02A = 0

drug_users = data.frame(id = drug_data$vdp_clean_id, 
                        J01_antibiotics = drug_data[, c("J01")],
                        A07A_antibiotics = drug_data_l3[, "A07A"],
                        S01A_antibiotics = drug_data_l3[, "S01A"],
                        A07_antidiarrheals = drug_data[, "A07"],
                        A06_anticonstipation = drug_data[, "A06"] )
colnames(drug_users) = c("id","J01_antibiotics","A07A_antibiotics","S01A_antibiotics","A07_antidiarrheals","A06_anticonstipation")

#################################
## 10) Prune the meta data down
##      to the necessary features
##      
#################################
######################################
## make a smoking never ever current variable
######################################
smoking_NEC = clinical_data[, c("age_when_started_to_smoke",
                                "average_cigarettes_consumption_per_day_last3months",
                                "average_cigarillos_consumption_per_day_last3months",
                                "average_cigars_consumption_per_day_last3months",
                                "average_pipe_tobacco_consumption_per_day_last3months",
                                "did_you_stop", 
                                "do_you_inhale", 
                                "have_quit",
                                "have_smoked_foafull_year",
                                "have_smoked_last3month") ]
NEC = ifelse(smoking_NEC$age_when_started_to_smoke == "Niet van toepassing", "never", "ever")
table(smoking_NEC$have_smoked_last3month, smoking_NEC$age_when_started_to_smoke)
w = which(smoking_NEC$have_smoked_last3month == 1)
NEC[w] = "current"
clinical_data$smoking_NEC = NEC
clinical_data$smoking_NEC = factor(clinical_data$smoking_NEC, levels = c("never","ever","current"))
mapping = c("never" = 0, "ever" = 1, "current" = 2)
b = mapping[clinical_data$smoking_NEC]
clinical_data = clinical_data |> mutate(smoking_NEC.r = b, .after = smoking_NEC)

######################################
## average food consumption last week
######################################
w = grep("average_consumption", colnames(clinical_data))
acv = clinical_data[,w]
for(i in 1:ncol(acv)){
  w = which(acv[,i] == "Tussen 5 en 10"); if(length(i)>0){ acv[w,i] = "6" }
  w = which(acv[,i] == "Meer dan 10"); if(length(i)>0){ acv[w,i] = "7" }
  w = which(acv[,i] == "Niet van toepassing"); if(length(i)>0){ acv[w,i] = "0" }
  acv[,i] = as.numeric( unlist(acv[,i]) )
}
acv = acv |> mutate(vdp_ids = clinical_data$vdp_ids, .before = colnames(acv)[1])

## Features to keep
features_2_keep = c("vdp_ids", 
                    "fasting", "ethnicity",
                    "gender_mq","gender_online",
                    "age_mq","age_online", 
                    "length","height","weight_mq","weight_online",
                    "BMI_mq","BMI_online", 
                    "hip_circumference","waist_circumference",
                    "smoking_NEC","smoking_NEC.r", "house_hold_monthl_net_income",
                    
                    "bristol_stool_score", "previous_relief", "time_of_sampling", 
                    "date_of_sampling", "abdominal_cramps_last_week", "average_defecations_per_day",
                    "bloating_last_week", "constipated_last_week", "defecation_days_count", 
                    "diarrhea_last_week", "flatulence_last_week", "have_used_diarrhea_inhibitors_last_week",
                    "intestinal_rumbling_last_week", "non_empty_bowels_feeling_after_defecation_last_week",
                    "used_laxatives_last_week",
                    
                    "BP.sys", "BP.dia", "sleeping_hours_per_day","exercising_habits",
                    
                    # "Enterotype",
                    "exclude_subject","exclusion_gender_mismatch",
                    "exclusion_age_mismatch","exclusion_BMI_mismatch",
                    
                    "is_vegetarian","is_vegan", "is_followingadiet",
                    # "last_weeks_diets", "special_diets", "diet_followed_for",
                    "followaspecial_dietarty_habit", "followedadiet_last_week",
                    "lactose_free_diet", "gluten_free_diet",
                    "weight_loss_diet", "fiber_rich_diet", "limited_fat_diet",
                    
                    "have_diabetes","diabetes_age_diagnosed","diabetes_type",
                    "have_high_blood_pressure", "have_high_blood_cholesterol", 
                    "vigorous_activities", 
                    "physical_function_score",  "social_function_score",
                    "physical_role_limited_score",  "emotional_role_limited_score",
                    "mental_health_score", "vitality_score", "pain_score",
                    "general_health_score"
                    )

meta_data = clinical_data[, features_2_keep]

###############################
## Add ACV to meta_data
###############################
m = match(meta_data$vdp_ids, acv$vdp_ids)
meta_data = cbind(meta_data, acv[, -1])

###############################
## Add PCS to meta_data
###############################
m = match(meta_data$vdp_ids, rownames(pcs))
meta_data = cbind( meta_data, pcs[m,])

###############################
## Add GRS data to meta_data
###############################
m = match(meta_data$vdp_ids, rownames(new_grs_data))
meta_data = cbind( meta_data, new_grs_data[m,-1])

###############################
## Add drug data to meta_data
###############################
m = match(meta_data$vdp_ids, drug_users$id )
meta_data = cbind( meta_data, drug_users[m,-1])


######################################
##
## Edit some of the variables
##
######################################

###################
### Ethnicity
###################
meta_data$ethnicity = gsub("blank/Zuid-Europees, mediterraan of Arabisch", "SEuro_Medit_Arabian", meta_data$ethnicity)
meta_data$ethnicity = gsub("blank/Oost- en West-Europees", "E_W_European", meta_data$ethnicity)
meta_data$ethnicity = gsub("Aziatisch", "Asian", meta_data$ethnicity)
w = grep("belg", meta_data$ethnicity); meta_data$ethnicity[w] = "E_W_European"
meta_data$ethnicity = gsub("Duits", "E_W_European", meta_data$ethnicity)
meta_data$ethnicity = gsub("nl", "E_W_European", meta_data$ethnicity)
meta_data$ethnicity = gsub("NL", "E_W_European", meta_data$ethnicity)
meta_data$ethnicity = gsub("nihil", "nothing", meta_data$ethnicity)
meta_data$ethnicity = gsub("Xxxx", "nothing", meta_data$ethnicity)
meta_data$ethnicity = gsub("bunderstraat 115 1745 Opwijk", "nothing", meta_data$ethnicity)
meta_data$ethnicity = gsub("zwart/negro&iuml;de", "black", meta_data$ethnicity)
meta_data$ethnicity = gsub("Xxxx", "nothing", meta_data$ethnicity)
meta_data$ethnicity = gsub( "nothing", "decline", meta_data$ethnicity)

###################
### house_hold_monthl_net_income
###################
meta_data$house_hold_monthl_net_income = gsub("Ik weet het niet", "I_dont_know", meta_data$house_hold_monthl_net_income)
meta_data$house_hold_monthl_net_income = gsub("Ik wil hier liever geen antwoord op geven", "rather_not_answer", meta_data$house_hold_monthl_net_income)
meta_data$house_hold_monthl_net_income = gsub("Minder dan â‚¬ 750", "less_than_750", meta_data$house_hold_monthl_net_income)
meta_data$house_hold_monthl_net_income = gsub("â‚¬ 750 - â‚¬ 1000", "750_to_1000", meta_data$house_hold_monthl_net_income)
meta_data$house_hold_monthl_net_income = gsub("â‚¬ 1000 - â‚¬ 1500", "1000_to_1500", meta_data$house_hold_monthl_net_income)
meta_data$house_hold_monthl_net_income = gsub("â‚¬ 1500 - â‚¬ 2000", "1500_to_2000", meta_data$house_hold_monthl_net_income)
meta_data$house_hold_monthl_net_income = gsub("â‚¬ 2000 - â‚¬ 2500", "2000_to_2500", meta_data$house_hold_monthl_net_income)
meta_data$house_hold_monthl_net_income = gsub("â‚¬ 2500 - â‚¬ 3000", "2500_to_3000", meta_data$house_hold_monthl_net_income)
meta_data$house_hold_monthl_net_income = gsub("â‚¬ 3000 - â‚¬ 3500", "3000_to_3500", meta_data$house_hold_monthl_net_income)
meta_data$house_hold_monthl_net_income = gsub("Meer dan â‚¬ 3500", "more_than_3500", meta_data$house_hold_monthl_net_income)
meta_data$house_hold_monthl_net_income = factor(meta_data$house_hold_monthl_net_income, levels = c("I_dont_know", "rather_not_answer", "less_than_750", 
                                                                                                      "1000_to_1500", "1500_to_2000", "2000_to_2500", "2500_to_3000", 
                                                                                                      "3000_to_3500", "more_than_3500") )
#### Recode to a numeric
mapping = c("I_dont_know" = NA, "rather_not_answer" = NA, 
            "less_than_750" = 1, "750_to_1000" = 2,  "1000_to_1500" = 3, "1500_to_2000" = 4, 
            "2000_to_2500" = 5, "2500_to_3000" = 6, "3000_to_3500" = 7, "more_than_3500" = 8  )
house_hold_monthl_net_income.r = mapping[meta_data$house_hold_monthl_net_income] 

meta_data = meta_data |> mutate(house_hold_monthl_net_income.r = house_hold_monthl_net_income.r, .after = house_hold_monthl_net_income)

###################
### previous_relief
###################
table( meta_data$previous_relief )
meta_data$previous_relief = gsub("Minder dan 6 uur geleden", "less_than_6_hours", meta_data$previous_relief)
meta_data$previous_relief = gsub("Tussen 6 en 12 uur geleden", "6_to_12_hours", meta_data$previous_relief)
meta_data$previous_relief = gsub("Tussen 12 en 18 uur geleden", "12_to_18_hours", meta_data$previous_relief)
meta_data$previous_relief = gsub("Tussen 18 en 24 uur geleden", "18_to_24_hours", meta_data$previous_relief)
meta_data$previous_relief = gsub("Tussen 24 en 36 uur geleden", "24_to_36_hours", meta_data$previous_relief)
meta_data$previous_relief = gsub("Tussen 36 en 48 uur geleden", "36_to_48_hours", meta_data$previous_relief)
meta_data$previous_relief = gsub("Langer dan 48 uur geleden", "more_than_48_hours", meta_data$previous_relief)

mapping = c("less_than_6_hours" = 1, "6_to_12_hours" = 2, "12_to_18_hours" = 3, 
            "18_to_24_hours" = 4, "24_to_36_hours" = 5, "36_to_48_hours" = 6, 
            "more_than_48_hours" = 7  )
previous_relief.r = mapping[meta_data$previous_relief] 

meta_data = meta_data |> mutate(previous_relief.r = previous_relief.r, .after = previous_relief)

###################
### time_of_sampling
###################
meta_data$time_of_sampling = gsub("9uur","09:00",meta_data$time_of_sampling)
meta_data$time_of_sampling = gsub("\\.",":",meta_data$time_of_sampling)
meta_data$time_of_sampling = lubridate::hm(meta_data$time_of_sampling)

###################
### date_of_sampling
###################
meta_data$date_of_sampling = lubridate::date( lubridate::ymd_hms(meta_data$date_of_sampling, tz = "CET") )

###################
## average_defecations_per_day
###################
x = meta_data$average_defecations_per_day
x = gsub("Meer dan 5","more_than_5", x)
meta_data$average_defecations_per_day = x
mapping = c("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5,"more_than_5" = 6)
b = mapping[x]
meta_data = meta_data |> mutate(average_defecations_per_day.r = b, .after = average_defecations_per_day)

############################
## A function to help speed 
## up uniform conversions
############################
flemish_to_english = function(x){
  x = gsub("dagen","days", x)
  x = gsub("dag","day", x)
  x = gsub("dag","day", x)
  x = gsub("Ik had hier de afgelopen week geen last van","more_than_7_days_ago", x)
  x = gsub("Ik heb hier zelden of nooit last van","not_an_issue", x)
  x = gsub("Ik gebruik deze zelden of nooit", "not_an_issue",x)
  x = gsub("Ik gebruikte er de afgelopen week geen", "more_than_7_days_ago",x)
  x = gsub(" ","_", x)
  return(x)
}
############################
## Uniform mapping
############################
mapping = c("1_day" = 9, "2_days" = 8, "3_days" = 7, 
            "4_days" = 6, "5_days" = 5, "6_days" = 4, 
            "7_days" = 3 , "more_than_7_days_ago" = 2, "not_an_issue" = 1 )

###################
### abdominal_cramps_last_week
###################
a = flemish_to_english(meta_data$abdominal_cramps_last_week)
meta_data$abdominal_cramps_last_week = a
b = mapping[a] 
meta_data = meta_data |> mutate(abdominal_cramps_last_week.r = b, .after = abdominal_cramps_last_week)


###################
### bloating_last_week
###################
a = flemish_to_english(meta_data$bloating_last_week)
meta_data$bloating_last_week = a
b = mapping[a] 
meta_data = meta_data |> mutate(bloating_last_week.r = b, .after = bloating_last_week)

###################
### constipated_last_week
###################
a = flemish_to_english(meta_data$constipated_last_week)
meta_data$constipated_last_week = a
b = mapping[a] 
meta_data = meta_data |> mutate(constipated_last_week.r = b, .after = constipated_last_week)

###################
### diarrhea_last_week
###################
a = flemish_to_english(meta_data$diarrhea_last_week)
meta_data$diarrhea_last_week = a
b = mapping[a] 
meta_data = meta_data |> mutate(diarrhea_last_week.r = b, .after = diarrhea_last_week)

###################
### flatulence_last_week
###################
a = flemish_to_english(meta_data$flatulence_last_week)
meta_data$flatulence_last_week = a
b = mapping[a] 
meta_data = meta_data |> mutate(flatulence_last_week.r = b, .after = flatulence_last_week)

###################
### have_used_diarrhea_inhibitors_last_week
###################
a = flemish_to_english(meta_data$have_used_diarrhea_inhibitors_last_week)
meta_data$have_used_diarrhea_inhibitors_last_week = a
b = mapping[a] 
meta_data = meta_data |> mutate(have_used_diarrhea_inhibitors_last_week.r = b, .after = have_used_diarrhea_inhibitors_last_week)

###################
### intestinal_rumbling_last_week
###################
a = flemish_to_english(meta_data$intestinal_rumbling_last_week)
meta_data$intestinal_rumbling_last_week = a
b = mapping[a] 
meta_data = meta_data |> mutate(intestinal_rumbling_last_week.r = b, .after = intestinal_rumbling_last_week)

###################
### non_empty_bowels_feeling_after_defecation_last_week
###################
a = flemish_to_english(meta_data$non_empty_bowels_feeling_after_defecation_last_week)
meta_data$non_empty_bowels_feeling_after_defecation_last_week = a
b = mapping[a] 
meta_data = meta_data |> mutate(non_empty_bowels_feeling_after_defecation_last_week.r = b, .after = non_empty_bowels_feeling_after_defecation_last_week)

###################
### used_laxatives_last_week
###################
a = flemish_to_english(meta_data$used_laxatives_last_week)
meta_data$used_laxatives_last_week = a
b = mapping[a] 
meta_data = meta_data |> mutate(used_laxatives_last_week.r = b, .after = used_laxatives_last_week)

###################
### exercising_habits
###################
x = meta_data$exercising_habits
x = gsub("ik doe weinig of nooit aan sport","rarely_never", x)
x = gsub("ik-sport 1-3 keer per week","1_to_3_times_a_week", x)
x = gsub("ik sport 3-5 keer per week","3_to_5_times_a_week", x)
x = gsub("ik sport 6-7 keer per week","6_to_7_times_a_week", x)
x = gsub("sport dagelijks/doe zwaar lichamelijk werk","daily_or_physical_work", x)
x = gsub("ik ","", x)
meta_data$exercising_habits = x

mapping = c("rarely_never" = 0, "1_to_3_times_a_week" = 1, 
            "3_to_5_times_a_week" = 2, "6_to_7_times_a_week" = 3, 
            "daily_or_physical_work" = 4)
b = mapping[x] 
meta_data = meta_data |> mutate(exercising_habits.r = b, .after = exercising_habits)

###################
### followaspecial_dietarty_habit
###################
x = meta_data$followaspecial_dietarty_habit
x = gsub("Nee","no", x)
x = gsub("Ja, flexibel","yes_flexible", x)
x = gsub("Ja, strikt","yes_strict", x)
meta_data$followaspecial_dietarty_habit = x
mapping = c("no" = 0, 
            "yes_flexible" = 1, 
            "yes_strict" = 2)
b = mapping[x] 
meta_data = meta_data |> mutate(followaspecial_dietarty_habit.r = b, .after = followaspecial_dietarty_habit)


###################
### followedadiet_last_week
###################
x = meta_data$followedadiet_last_week
x = gsub("Nee","no", x)
x = gsub("Ja, enkele dagen","yes_few_days", x)
x = gsub("Ja, de hele week","yes_all_week", x)
meta_data$followedadiet_last_week = x
mapping = c("no" = 0, 
            "yes_few_days" = 1, 
            "yes_all_week" = 2)
b = mapping[x] 
meta_data$followedadiet_last_week = x
meta_data = meta_data |> mutate(followedadiet_last_week.r = b, .after = followedadiet_last_week)


###################
### have_diabetes
###################
x = meta_data$have_diabetes
x = gsub("Neen","no", x)
x = gsub("Ja","yes", x)
meta_data$have_diabetes = x
mapping = c("no" = 0, 
            "yes" = 1 )
b = mapping[x] 
meta_data = meta_data |> mutate(have_diabetes.r = b, .after = have_diabetes)

###################
### diabetes_type
###################
x = meta_data$diabetes_type
x = gsub("type 1 \\(jeugddiabetes, meestal van kinds af aan\\)","typeI", x)
x = gsub("type 2 \\(ouderdomsdiabetes, meestal op latere leeftijd ontstaan\\)","typeII", x)
x[which( !x %in% c("typeI", "typeII", NA)  )] = "other"

meta_data$diabetes_type = x


###################
### have_high_blood_pressure
###################
x = meta_data$have_high_blood_pressure
x = gsub("Neen","no", x)
x = gsub("Ja","yes", x)
meta_data$have_high_blood_pressure = x
mapping = c("no" = 0, 
            "yes" = 1 )
b = mapping[x] 
meta_data = meta_data |> mutate(have_high_blood_pressure.r = b, .after = have_high_blood_pressure)

###################
### have_high_blood_cholesterol
###################
x = meta_data$have_high_blood_cholesterol
x = gsub("Neen","no", x)
x = gsub("Ja","yes", x)
meta_data$have_high_blood_cholesterol = x
mapping = c("no" = 0, 
            "yes" = 1 )
b = mapping[x] 
meta_data = meta_data |> mutate(have_high_blood_cholesterol.r = b, .after = have_high_blood_cholesterol)

###################
### vigorous_activities
###################
x = meta_data$vigorous_activities
x = gsub("Nee, helemaal niet beperkt","no_limitations", x)
x = gsub("Ja, een beetje beperkt","yes_limited", x)
x = gsub("Ja, ernstig beperkt","yes_severly_limited", x)
meta_data$vigorous_activities = x
mapping = c("no_limitations" = 2, 
            "yes_limited" = 1,
            "yes_severly_limited" = 0)
b = mapping[x] 
meta_data = meta_data |> mutate(vigorous_activities.r = b, .after = vigorous_activities)


#######################################
## 11) Build the Enterotype binary traits
##
#######################################
a = mt_data$EnterotypeClass
a[a == "Bact1"] = "1"
a[a != "1"] = "0"
mt_data = mt_data |> mutate(Enterotype_Bact1 = a, .before = MDS1)
###
a = mt_data$EnterotypeClass
a[a == "Bact2"] = "1"
a[a != "1"] = "0"
mt_data = mt_data |> mutate(Enterotype_Bact2 = a, .before = MDS1)
##
a = mt_data$EnterotypeClass
a[a == "Prev"] = "1"
a[a != "1"] = "0"
mt_data = mt_data |> mutate(Enterotype_Prev = a, .before = MDS1)
##
a = mt_data$EnterotypeClass
a[a == "Rum"] = "1"
a[a != "1"] = "0"
mt_data = mt_data |> mutate(Enterotype_Rum = a, .before = MDS1)

#################################
## 9) Define WHRadjBMI
#################################
cat(paste0("3. Defining WHR and WHRadjBMI\n"))
# alt_clinical_data$whr = alt_clinical_data$Buikoptrek/ alt_clinical_data$Heupomtrek
##
meta_data = meta_data |> mutate(WHR = waist_circumference / hip_circumference, .before = "smoking_NEC" )

fit = lm(WHR ~ BMI_online, data = meta_data)
res = fit$residuals
m = match(rownames(meta_data), names(res))
meta_data = meta_data |> mutate(WHRadjBMI = res[m], .before = "smoking_NEC" )


#################################
## 10) remove samples with no MT data
#################################
cat(paste0("4. Removing samples with NO microbiome data\n"))
missing = apply( mt_data[, c(3:509)] , 1, function(x){ sum( is.na(x) ) })
w = which(missing>0)  ## n = 2
##
if(length(w)>0){
  meta_data = meta_data[-w, ]
  mt_data = mt_data[-w, ]
  metabolon_data = metabolon_data[-w, ]
  batch_data = batch_data[-w,]
  ### Other Source data frames
  clinical_data = clinical_data[-w, ]
  gwased_mt_data = gwased_mt_data[-w, ]
  grs_data = grs_data[-w, ]
  new_grs_data = new_grs_data[-w, ]
  blood_traits = blood_traits[-w, ]
  
}

#################################
## 11) Pull the data together
#################################
cat(paste0("5. Pull everything together\n"))
## PULL EVERYTHING TOGETEHR
study_data = list()
study_data$meta_data = meta_data
study_data$micro_data = mt_data
study_data$batch_data = batch_data
study_data$metabolon_data = metabolon_data
study_data$metabolon_feature_data = metabolon_feature_data
### Other SOURCE Data
study_data$clinical_data = clinical_data
study_data$blood_traits = blood_traits
study_data$gwased_mt_data = gwased_mt_data
study_data$grs_data = grs_data
study_data$new_grs_data = new_grs_data

#################################
## 12) SAVE data to Rdata file
#################################
cat(paste0("6. write to file\n"))
n = paste0( project_data_dir , "study_data_20240911_v1.Rdata")
# save(study_data, file = n)


#################################
##
## 13) IMPUTATION
##
#################################
cat(paste0("7. Imputation\n"))
## how much missingness is there?
vars = c("age_online", "age_mq", 
         "gender_online", "gender_mq", 
         "weight_online", "weight_mq",
         "height",
         "BMI_online","BMI_mq",
         "hip_circumference","waist_circumference", "WHR", "WHRadjBMI")
         
missing_count = apply( meta_data[, vars ], 2, function(x){sum(is.na(x))} )
# age_online          age_mq       gender_online           gender_mq       weight_online 
# 56                  29                   8                  30                  19 
# weight_mq          height          BMI_online              BMI_mq   hip_circumference 
# 56                  26                  26                  60                 144 
# waist_circumference                 WHR           WHRadjBMI 
# 145                 154                 177 
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
dim(working_mt) ## 833

## remove AB and PA
w = grep("PA", colnames(working_mt))
working_mt = working_mt[,-w]
dim(working_mt) ## 619

w = grep("AB", colnames(working_mt))
working_mt = working_mt[,-w]
dim(working_mt) ## 509

## Identify microbioal traits with the least zero
zero_prop = apply(working_mt, 2, function(x) { sum(x == 0, na.rm = TRUE)/length( na.omit(x) )   } )
w = which(zero_prop < 0.1 )
working_mt = working_mt[,w]
dim(working_mt) ## 413

## keep only those with a mean of 50
mean_est = apply(working_mt, 2, function(x) {  mean( as.numeric(x), na.rm = TRUE)   } )
w = which(mean_est > 50)
working_mt = working_mt[,w]
for(i in 1:ncol(working_mt)){ working_mt[,i] = as.numeric(working_mt[,i])}
dim(working_mt) ## 82

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
dim(working_mt) ## 48

## Evaluate tree cut
Cmat = cor(working_mt, method = "sp")
Dmat = as.dist(1 - abs(Cmat) )
nj = hclust(Dmat, method = "complete")
# plot(nj, hang = -1)

#######################
## Blood traits
#######################
mis = apply(blood_traits, 2, function(x){ sum(is.na(x)) } )
w = names( which(mis<= min(mis) ) )
working_blood_traits = blood_traits[,  w]

###########################
### IMPUTATION DATA SET
###########################
anthro_vars = c("age_online", "age_mq", 
                "gender_online", "gender_mq", 
                "weight_online", "weight_mq",
                "height",
                "BMI_online","BMI_mq",
                "hip_circumference","waist_circumference", 
                "WHR",
                "WHRadjBMI")

############
imputation_data = cbind(meta_data[, anthro_vars],
                        working_blood_traits, 
                        working_mt, 
                        working_metabolon)
############
variables_2_test_for_imputation = colnames(imputation_data)[-c(1:13)]

############################
##
## *** Perform Imputation ***
##
###########################
vars = colnames(imputation_data)[1:13]


########################
### AGE IMPUTATION
########################
cat(paste0("7a. AGE Imputation\n"))
imputation = 0(wdata = imputation_data, 
                            var_2_impute = "age_mq",
                            r2_threshold = 0.05,
                            variables_2_test_for_imputation = variables_2_test_for_imputation )
## sumstats
unlist( imputation[2:4] )
## average imputation error
imputation[5] ##  6.14 (95CI 0.558-21.664)


## Add Imputation to data set
meta_data = meta_data |> mutate(imputed_age = meta_data$age_mq, .before = "length")
m = match( names(imputation[[1]]) , rownames(meta_data) )
meta_data$imputed_age[m] = as.character( imputation[[1]] )


########################
### SEX IMPUTATION
########################
cat(paste0("7b. SEX Imputation\n"))
sex_imputation = knn_impute(wdata = imputation_data, 
           var_2_impute = "gender_online",
           r2_threshold = 0.05,
           variables_2_test_for_imputation = variables_2_test_for_imputation )
## sumstats
unlist( sex_imputation[2:4] )
## average imputation error
sex_imputation[5] ## 0% (95CI 0-0.25)
 
## Add Imputation to data set
meta_data = meta_data |> mutate(imputed_sex = meta_data$gender_online, .before = "age_mq")
m = match( names(sex_imputation[[1]]) , rownames(meta_data) )
meta_data$imputed_sex[m] = as.character( sex_imputation[[1]] )

########################
### Weight IMPUTATION
########################
cat(paste0("7c. WEIGHT Imputation\n"))
imputation = knn_impute(wdata = imputation_data, 
                            var_2_impute = "weight_online",
                            r2_threshold = 0.05,
                            variables_2_test_for_imputation = variables_2_test_for_imputation )
## sumstats
unlist( imputation[2:4] )
## average imputation error
imputation[5] ## 7.62 (95% CI 0.998-25.53)

## Add Imputation to data set
meta_data = meta_data |> mutate(imputed_weight_online = meta_data$weight_online, .before = "BMI_mq")
m = match( names(imputation[[1]]) , rownames(meta_data) )
meta_data$imputed_weight_online[m] = as.numeric( as.character( imputation[[1]] ) )



########################
### Height IMPUTATION
########################
cat(paste0("7d. HEIGHT Imputation\n"))
## remove the blood cell traits as one individual who is missing height has no blood cell trait data
w = which( variables_2_test_for_imputation %in% c("WHRadjBMI",colnames(blood_traits) ) )

imputation = knn_impute(wdata = imputation_data, 
                        var_2_impute = "height",
                        r2_threshold = 0.05,
                        variables_2_test_for_imputation = variables_2_test_for_imputation[-w] )
## sumstats
unlist( imputation[2:4] )
## average imputation error
imputation[5] ## 5.11 (95% CI 0.28-16.22)

## Add Imputation to data set
meta_data = meta_data |> mutate(imputed_height = meta_data$height, .before = "weight_mq")
m = match( names(imputation[[1]]) , rownames(meta_data) )
meta_data$imputed_height[m] = as.numeric( as.character( imputation[[1]] ) )


########################
### BMI IMPUTATION
########################
cat(paste0("7e. BMI Imputation\n"))
## Remove blood cell traits from imputation variables
w = which( variables_2_test_for_imputation %in% c("WHRadjBMI",colnames(blood_traits) ) )

imputation = knn_impute(wdata = imputation_data, 
                        var_2_impute = "BMI_online",
                        r2_threshold = 0.05,
                        variables_2_test_for_imputation = variables_2_test_for_imputation[-w] )
## sumstats
unlist( imputation[2:4] )
## average imputation error
imputation[5] ## 2.58 (95% CI 0.259-9.36)

## Add Imputation to data set
meta_data = meta_data |> mutate(imputed_BMI_online = meta_data$BMI_online, .before = "hip_circumference")
m = match( names(imputation[[1]]) , rownames(meta_data) )
meta_data$imputed_BMI_online[m] = round(as.numeric( as.character( imputation[[1]] ) ), d = 5)



########################
### hip_circumference IMPUTATION
########################
cat(paste0("7f. HIP Cir. Imputation\n"))
## Remove blood cell traits from imputation variables
w = which( variables_2_test_for_imputation %in% c("WHRadjBMI",colnames(blood_traits) ) )

imputation = knn_impute(wdata = imputation_data, 
                        var_2_impute = "hip_circumference",
                        r2_threshold = 0.05,
                        variables_2_test_for_imputation = variables_2_test_for_imputation[-w] )
## sumstats
unlist( imputation[2:4] ) 
## average imputation error
imputation[5] ## 6.133 (95% CI 0.155-23.93)

## Add Imputation to data set
meta_data = meta_data |> mutate(imputed_hip_circumference = meta_data$hip_circumference, .before = "waist_circumference")
m = match( names(imputation[[1]]) , rownames(meta_data) )
meta_data$imputed_hip_circumference[m] =  as.numeric( as.character( imputation[[1]] ) )


########################
### waist_circumference IMPUTATION
########################
cat(paste0("7g. WAIST Cir. Imputation\n"))
w = which( variables_2_test_for_imputation %in% c("WHRadjBMI",colnames(blood_traits) ) )

imputation = knn_impute(wdata = imputation_data, 
                        var_2_impute = "waist_circumference",
                        r2_threshold = 0.05,
                        variables_2_test_for_imputation = variables_2_test_for_imputation[-w] )
## sumstats
unlist( imputation[2:4] ) 
## average imputation error
imputation[5] ## 7.48 (95% CI 0.267-26.87)

## Add Imputation to data set
meta_data = meta_data |> mutate(imputed_waist_circumference = meta_data$waist_circumference, .before = "WHR")
m = match( names(imputation[[1]]) , rownames(meta_data) )
meta_data$imputed_waist_circumference[m] =  as.numeric( as.character( imputation[[1]] ) )


########################
### WHR IMPUTATION
########################
cat(paste0("7h. WHR Imputation\n"))
w = which( variables_2_test_for_imputation %in% c("WHRadjBMI",colnames(blood_traits) ) )

imputation = knn_impute(wdata = imputation_data, 
                        var_2_impute = "WHR",
                        r2_threshold = 0.025,
                        variables_2_test_for_imputation = variables_2_test_for_imputation[-w] )
## sumstats
unlist( imputation[2:4] ) ## of vars used = 93
## average imputation error
imputation[5] ## 0.0638 (95% CI 0.0031-0.259)

## Add Imputation to data set
meta_data = meta_data |> mutate(imputed_WHR = meta_data$WHR, .before = "WHRadjBMI")
m = match( names(imputation[[1]]) , rownames(meta_data) )
meta_data$imputed_WHR[m] =  as.numeric( as.character( imputation[[1]] ) )


########################
### WHRadjBMI IMPUTATION
########################
cat(paste0("8. Est WHRadjBMI\n"))
fit = lm(imputed_WHR ~ imputed_BMI_online, data = meta_data)
res = fit$residuals
m = match(rownames(meta_data), names(res))
meta_data = meta_data |> mutate(imputed_WHRadjBMI = res[m], .before = "smoking_NEC")


#################################
## 14) Pull the data together
#################################
cat(paste0("9. Pull Data Together\n"))

## Redefine Study meta_data
study_data$meta_data = meta_data

#################################
## 15) SAVE data to Rdata file
#################################
cat(paste0("10. Write to file\n"))
n = paste0( project_data_dir , "study_data_20240911_v1.1.Rdata")
save(study_data, file = n)


