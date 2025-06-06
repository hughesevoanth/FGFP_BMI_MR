################################################################
## FGFP MR: BMI-->Microbiome
## -- Load data and functions for scripts --
## by: David Hughes
## date: March 16th 2022
################################################################

############################################
## Load parameters from parameter file
############################################
source("parameters/pfile.sh")


############################################
## Load functions
############################################
source( "functions/taxa_2_keep.R" )
source( "functions/build_mts.R" )
source( "functions/rntransform.R" )
source( "functions/obsmr.R" )
source( "functions/est_beta_varexp.R" )
source( "functions/id_outliers.R" )
source( "functions/ivregfit.R" )
source( "functions/lmfit.R" )
source( "functions/lmfit_v3.R" )
source( "functions/normW.R" )
source( "functions/obsmr_2_long.R" )
source( "functions/dh_forrest_plot.R" )
source( "functions/ivtoolsfit.R" )
source( "functions/ivtoolsfit_2_long.R" )
source( "functions/ztest.R" )


############################################
## load the source data
############################################
f = paste0(project_data_dir, "study_data_20240917_v1.2.Rdata")
load( f )

