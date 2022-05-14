ivtoolsfit = function( wdata,
                       outcome,
                       exposure,
                       instrument,
                       covariates = NA,
                       outcome_model_family = "gaussian",
                       exposure_model_family = "gaussian",
                       weights = NA,
                       rnt_outcome = FALSE){

  ### Define the model data
  model_variables = na.omit( c(outcome, exposure, instrument, covariates, weights) )
  mod_data = na.omit( wdata[, c(model_variables)] )
  
  ############################################
  ## I. rank normalize the dependent|outcome ?
  ############################################
  if(outcome_model_family == "gaussian" & rnt_outcome == TRUE){
    set.seed(2022)
    mod_data[, outcome] = rntransform( mod_data[, outcome] )
  }
  
  ################################################
  ## II. IF OUTCOME binomial make sure outcome is a factor
  ################################################
  if(outcome_model_family == "binomial"){
    mod_data[, outcome] = as.factor(mod_data[, outcome])
  }

  ########################
  ### III. Build STAGE ONE
  ##      model
  ########################
  if( length(na.omit(covariates))>0 ){
    S1_form = formula(paste0( exposure, " ~ ", paste0( covariates, collapse = " + ") , "+", instrument ) )
  } else {
    S1_form = formula( paste0( exposure, " ~ ",  instrument ) )
  }

  #########################
  ## IV. RUN STAGE ONE
  #########################
  ## STAGE ONE
  if( is.na( weights) ){
    fitX.LZ <- glm(formula = S1_form , data = mod_data, family = exposure_model_family)
  } else {
    w = grep( weights, colnames(mod_data) ); colnames(mod_data)[w] = "study_weights"
    fitX.LZ <- glm(formula = S1_form , data = mod_data, weights = study_weights, family = exposure_model_family)
  }
  
  #########################
  ## IVa.
  ## Extract Exposure on Instrument (ei)
  ## Summary Statistics
  #########################
  res = residuals(fitX.LZ)
  ei_n = length(res)
  names(ei_n) = "ei_n"
  
  ei_summary = summary(fitX.LZ)
  ei_coef = ei_summary$coefficients[instrument, ]
  ei_stat_n = paste0("ei_", gsub(" ","",names(ei_coef)[3] ) )
  names( ei_coef ) = c("ei_beta","ei_se",ei_stat_n,"ei_P")
  ## R2
  a = car::Anova(fitX.LZ, type = "II", test = "F")
  ei_Fstat = a[instrument,3]
  names(ei_Fstat) = "ei_Fstat"
  ss = a[,1]
  names(ss) = rownames(a)
  ei_eta_sq = ss[instrument] / sum(ss)
  names(ei_eta_sq) = "ei_eta_sq"
  ei_stats = c(ei_n, ei_coef, ei_Fstat, ei_eta_sq)
  
  ########################
  ### V Build STAGE TWO
  ##      model
  ########################
  if( length(na.omit(covariates))>0 ){
    S2_form = formula(paste0( outcome, " ~ ", paste0( covariates, collapse = " + ") ,"+" , exposure ) )
  } else {
    S2_form = formula( paste0( outcome, " ~ ",  exposure ) )
  }
  
  #########################
  ## VI. RUN STAGE TWO
  #########################
  ## STAGE ONE
  if( is.na( weights) ){
    fitY.LX <- glm(formula = S2_form , data = mod_data, family = outcome_model_family)
  } else {
    w = grep( weights, colnames(mod_data) ); colnames(mod_data)[w] = "study_weights"
    fitY.LX <- glm(formula = S2_form , data = mod_data, weights = study_weights, family = outcome_model_family)
  }
  
  #########################
  ## Va.
  ## Extract Observational 
  ## regress outcome on exposure (oe)
  ## Summary Statistics
  #########################
  ## Sample Size
  res = residuals(fitY.LX)
  oe_n = length(res)
  names(oe_n) = "oe_n"
  
  ## normality of residuals
  if(outcome_model_family == "gaussian"){
    oe_Wstat = normW(res)
    names(oe_Wstat) = "oe_W"
  }
  
  if(outcome_model_family == "binomial"){
    oe_Wstat = NA
    names(oe_Wstat) = "oe_W"
  }
  
  ## Breusch-Pagan test of homoskedasticity
  oe_Breusch_Pagan_P = lmtest::bptest(fitY.LX)$p.value
  names(oe_Breusch_Pagan_P) = "oe_Breusch_Pagan_P"
  
  ## Model summary Stats
  oe_summary = summary(fitY.LX)
  
  ## Model coefficients
  oe_coef = oe_summary$coefficients[exposure, ]
  oe_stat_n = paste0("oe_", gsub(" ","",names(oe_coef)[3] ) )
  names( oe_coef ) = c("oe_beta","oe_se",oe_stat_n,"oe_P")
  
  ## Type II ANOVA
  a = car::Anova(fitY.LX, type = "II", test = "F")
  
  ## Exposure Fstat
  oe_Fstat = a[exposure,3]
  names(oe_Fstat) = "oe_Fstat"
  
  ## WALD TEST
  # car::linearHypothesis(ivmod, paste0(exposure, "=0"), test = "F" )
  oe_Wald_F_test = lmtest::waldtest(fitY.LX, exposure, test = "F" )
  oe_Wald_F_test = c(oe_Wald_F_test$F[2], oe_Wald_F_test$`Pr(>F)`[2])
  names(oe_Wald_F_test) = c("oe_Wald_F", "oe_Wald_P")
  
  
  ## Variance Explained (eta-sq  or R2)
  ss = a[,1]
  names(ss) = rownames(a)
  
  ## Model Eta-sq
  w = which(names(ss) == "Residuals")
  oe_model_eta_sq = sum( ss[-w] ) / sum( ss )
  names(oe_model_eta_sq) = "oe_model_eta_sq"
  
  ## Exposure Eta-sq
  oe_exposure_eta_sq = ss[exposure] / sum(ss)
  names(oe_exposure_eta_sq) = "oe_exposure_eta_sq"
  
  
  oe_stats = c(oe_n, oe_coef, oe_Fstat, oe_Wald_F_test, oe_model_eta_sq, oe_exposure_eta_sq, oe_Breusch_Pagan_P)
  
  
  #########################
  ## VII. RUN MR model
  #########################
  fitIV_ts <- ivtools::ivglm( estmethod = "ts", fitX.LZ = fitX.LZ, fitY.LX = fitY.LX, data = mod_data, ctrl = FALSE  )
  
  ######################
  ### IV. Summary Stats
  ######################
  ## IV model summary
  s = summary(fitIV_ts)
  
  ## model coefficients
  coef = s$coefficients
  
  ## Report beta, se, t-value, and P-value
  MR_coef  = coef[exposure,]
  cn = gsub(" ","", names(MR_coef)[3] )
  names(MR_coef) = c("beta","se", cn ,"P")
  names(MR_coef) = paste0("MR_", names(MR_coef) )
    
  ## The re-ffitted IV model (Y.LX)
  ivmod = fitIV_ts$fitY.LX
  
  ## sample size in model
  res = residuals( ivmod )
  iv_n = length(res); names(iv_n) = "MR_n"
  
  ## normality of residuals
  if(outcome_model_family == "gaussian"){
    Wstat = normW(res)
    names(Wstat) = "MR_res_Wstat"
  } else {
    Wstat = NA
    names(Wstat) = "MR_res_Wstat"
  }
  
  ## Breusch-Pagan test of homoskedasticity
  Breusch_Pagan_P = lmtest::bptest(ivmod)$p.value
  names(Breusch_Pagan_P) = "MR_Breusch_Pagan_P"
  
  ## TYPE II ANOVA
  a = car::Anova(fitIV_ts$fitY.LX, type = "II", test.statistic = "F")
  
  ## Exposure Fstat
  MR_Fstat = a[exposure,3]
  names(MR_Fstat) = "MR_Fstat"
  
  ## WALD TEST
  # car::linearHypothesis(ivmod, paste0(exposure, "=0"), test = "F" )
  MR_Wald_F_test = lmtest::waldtest(ivmod, exposure, test = "F" )
  MR_Wald_F_test = c(MR_Wald_F_test$F[2], MR_Wald_F_test$`Pr(>F)`[2])
  names(MR_Wald_F_test) = c("MR_Wald_F", "MR_Wald_P")
  
  ## Variance Explained (eta-sq  or R2)
  sumsq = a[,1]; names(sumsq) = rownames(a)
  
  ## Model Eta-sq
  w = which(names(sumsq) == "Residuals")
  model_r2 = sum( sumsq[-w] ) / sum( sumsq )
  names(model_r2) = "MR_model_Rsq"
  
  ##  Genotype Predicted Exposure eta-sq
  dhat_r2 = sumsq[ exposure ]  / sum( sumsq )
  names(dhat_r2) = "MR_dhat_r2"
  
  MR_stats = c(iv_n, Wstat, MR_coef, model_r2, dhat_r2, MR_Fstat, MR_Wald_F_test, Breusch_Pagan_P  )
  
  ############################
  ## Data OUT
  ############################
  names(exposure) = "MR_exposure"
  names(outcome) = "MR_outcome"
  names(instrument) = "MR_instrument"
  ###
  out = c(outcome, exposure, instrument, 
          ei_stats,
          oe_stats,
          MR_stats )
  
  return(out)

}
