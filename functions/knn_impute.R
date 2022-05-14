################################################################
## A Function to impute the missing anthropomorphic traits
##  where possible
## by: David Hughes
## date: March 20th 2022
################################################################
knn_impute = function(wdata, ## data frame of data with all working data 2 use in imputation
                      var_2_impute, ## column name of variable to use in imputation
                      r2_threshold = 0.1,  ## R2 value to use as a minimum threshold for variables to use in imputation
                      variables_2_test_for_imputation ){
  
  ###############################################################
  ## Once I remove the samples with missingness at my variable of
  ## interest, what other samples have NO missing data?
  ## It is these samples with no other missingness that will be
  ## carried forward
  ###############################################################
  ## The samples with missing data at variable 2 impute
  w =  which( is.na( wdata[, var_2_impute]) )
  IDs_with_missing_var = rownames(wdata)[w]
  
  ## what other samples have missing data elsewhere in the data set?
  mis = apply(wdata[-w,], 1, function(x){ sum(is.na(x)) } )
  sample_2_keep = names( which(mis == 0) )
  
  count_of_samples_used_4_imputations = length(sample_2_keep)
  
  ### Redefine wdata
  m = match(c(IDs_with_missing_var, sample_2_keep), rownames(wdata))
  wdata = wdata[m, ]
  
  ## How many of my variable of interest are missing ?
  missing_var_count = sum( is.na(wdata[, var_2_impute]) )
  
  ### starting data set
  simdata = na.omit( wdata[, c(var_2_impute, variables_2_test_for_imputation)] )
  
  ## RNT the wdata
  simdata[,-1] = apply(simdata[,-1], 2, rntransform)
  
  ### identify top associated variables
  cors = t( sapply( colnames(simdata)[-1], function(v){ 
    form = formula(paste0(v, "~ ", var_2_impute))
    fit = lm(form, data = simdata)
    s = summary(fit)
    P = s$coef[2,4]
    R2 = s$r.squared
    out = c(v, R2, P)
    names(out) = c("var","r2","p")
    return(out)
  }) )
  cors = as.data.frame(cors)
  for(i in 2:3){cors[,i] = as.numeric(cors[,i]) }
  cors = cors %>% arrange(p)
  
  ## Identify those that explain >= 10% of the variation in the variable 2 impute
  w = which(cors$r2 >= r2_threshold)
  vars_4_imputations = cors$var[w]
  number_of_variables_used_in_imputation = length(vars_4_imputations)
  
  ############################
  ##
  ## Simulate Imputation
  ##
  ###########################
  binary = ifelse( length( na.omit( unique(wdata[, var_2_impute]) ) ) == 2, TRUE, FALSE )
  
  ## iterate over 100 samples
  var_sims = sapply(1:100, function(i){
    ## randomly sample individuals for prediction
    s = sort( sample(1:nrow(simdata), missing_var_count) )
    ## training data
    train = simdata[-s, vars_4_imputations]
    train_var = simdata[-s, var_2_impute ]
    ## test data
    test = simdata[s, vars_4_imputations ]
    test_var = simdata[s, var_2_impute ]
    ##  knn IMPUTATION
    impute_var <- knn(train=train, test=test, cl=train_var, k=round( sqrt(nrow(train)) ) )
    if(binary == TRUE){
      ## error rate in binary prediction
      out = sum( test_var != impute_var ) / length(impute_var)
    } else{
      ## error rate in quantita prediction
      delta = abs( test_var - as.numeric(as.character(impute_var)) )
      out = quantile(delta, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))  
    }
    return(out)
  })
  ###
  ###### REPORT ERROR ESTIMATES
  if(binary == TRUE){
    percentile_errors = quantile(var_sims, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
  } else {
    percentile_errors = apply( t(var_sims), 2, mean )
  }
  
  ############################
  ##
  ##  PERFORM Imputation
  ##
  ###########################
  ### starting data set
  impute_data = wdata[, c(var_2_impute, vars_4_imputations)] 
  
  ## RNT the wdata
  impute_data[,-1] = apply(impute_data[,-1], 2, rntransform)
  
  ## remove those individuals with missing var_2_impute data
  mis = which(is.na(impute_data[, var_2_impute]))
  w = which(colnames(impute_data) == var_2_impute)
  test = impute_data[ mis, -w ]
  
  ## For the samples that I want to impute the var_2_impute
  ## They must have COMPLETE DATA at the variables used for imputation.
  ## Which, if any, do not? Remove the from the test data frame.
  id_sam_with_mis = which( apply(test, 1, function(x){ sum(is.na(x))  } ) > 0 )
  number_of_individuals_I_could_not_impute = 0
  if(length(id_sam_with_mis) > 0){ 
    test = test[-id_sam_with_mis, ]  
    number_of_individuals_I_could_not_impute = length(id_sam_with_mis)
    }
  
  ## training data
  w = which(colnames(impute_data) == var_2_impute)
  train = impute_data[ -mis, -w ]
  train_pheno = impute_data[-mis, var_2_impute]
  
  ### perform the imputation
  Imputed_Var <- knn(train=train, test=test, cl=train_pheno, k=round( sqrt(nrow(train)) ) )
  names(Imputed_Var) = rownames(test)
  
  ## Data out
  out = list(imputed_Var = Imputed_Var ,
             number_of_samples_used_4_imputations = count_of_samples_used_4_imputations,
             number_of_variables_used_in_imputation = number_of_variables_used_in_imputation, 
             number_of_individuals_I_could_not_impute = number_of_individuals_I_could_not_impute,
             imputation_error = percentile_errors,
             CorMat = cors
  )
}
