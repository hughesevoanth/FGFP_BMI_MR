lmfit_v2 = function( wdata,
                  outcome = NA,
                  exposure = NA,
                  covariates = NA,
                  weights = NA,
                  rnt_outcome = FALSE,
                  typeIIanova = TRUE){

  ############################################
  ## I. define data set
  ############################################
  wdata = wdata[, na.omit(c(covariates, outcome, exposure, weights)) ]
  wdata = na.omit(wdata)
  
  ############################################
  ## I. rank normalize the outcome ?
  ############################################
  if(rnt_outcome == TRUE){
    wdata[, outcome] = rntransform( wdata[, outcome] )
  }

  ######################
  ### II. Define Linear model
  ######################
  if( length(na.omit(covariates)) > 0 ){
    form = formula(paste0(outcome, " ~ ", paste0( covariates, collapse = " + ") , " + ", exposure ))
  } else {
    form = formula(paste0(outcome, " ~ ", exposure ))
  }

  #########################
  ## III. RUN Generalized LINEAR MODEL
  #########################
  if( is.na(weights) ){
    lm_mod = lm(form, data = wdata )
  } else {
    ## re-name weights column
    w = grep( weights, colnames(wdata) ); colnames(wdata)[w] = "study_weights"
    lm_mod = lm(form, weights = study_weights, data = wdata)  
  }
  
  #########################
  ## IV. Summary Stats
  #########################
  ## sample size in model
  res = residuals(lm_mod)
  n = length(res)
  names(n) = "n"

  ## normality of residuals
  Wstat = normW(res)
  names(Wstat) = "W"
  
  ## Breusch-Pagan test of homoskedasticity
  Breusch_Pagan_P = lmtest::bptest(lm_mod)$p.value
  names(Breusch_Pagan_P) = "Breusch_Pagan_P"
  
  ## model summary
  s = summary(lm_mod)
  
  ## model coefficients
  lm_coef = s$coefficients
  
  ## Variance explained by model
  model_R2 = s$r.squared; names(model_R2) = "model_R2"
  model_adjR2 = s$adj.r.squared; names(model_adjR2) = "model_adjR2"
  
  ## Report beta, se, t-value, and P-value
  # lm_estimates = lm_coef[exposure, ]; names(lm_estimates) = c("beta","se","tval","P")
  lm_estimates = lm_coef[nrow(lm_coef), ]; names(lm_estimates) = c("beta","se","tval","P")
  
  ## extract sums of squares
  if(typeIIanova == FALSE){
    a = anova(lm_mod)
    ss = a[,2]
    names(ss) = rownames(a)
    eta_sq = ss/sum(ss); names(eta_sq) = paste0("etasq_", names(eta_sq))
    ## covariable P_values
    covar_P = a[-nrow(a),5]
    names(covar_P) = paste0( rownames(a)[-nrow(a)], "_P")
    
  } else {
    a = car::Anova(lm_mod, type = "II")
    ss = a[,1]
    names(ss) = rownames(a)
    eta_sq = ss/sum(ss, na.rm = TRUE); names(eta_sq) = paste0("etasq_", names(eta_sq))
    ## covariable P_values
    covar_P = a[-nrow(a),4]
    names(covar_P) = paste0( rownames(a)[-nrow(a)], "_P")
    
  }
  
  

  
  ## Sandwich co-variance matrix
  ## ** removes the assumption that residual errors have constant variance **
  # vcov_mat = sandwich::vcovHC(glm_mod, type = "HC") ## HC = White’s estimator
  vcov_mat = sandwich::vcovHC(lm_mod, type = "HC") ## HC = White’s estimator
  ## Sandwich SE's
  sandwich_se <- diag(vcov_mat)^0.5
  ## re-estimate p_values and make Coefficient table
  coef_table = lmtest::coeftest(lm_mod, vcov = vcov_mat)
  # sandwich_estimates = coef_table[exposure, ]
  sandwich_estimates = coef_table[nrow(coef_table), ]
  names(sandwich_estimates) = paste0("sw_", c("beta","se","zval","P"))
  
  ######################
  ## V. Linear model data out
  ######################
  names(exposure) = "exposure"
  names(outcome) = "outcome"
  
  # glm_out = c(exposure, outcome, n, Wstat, Breusch_Pagan_P, devexp_by_model, devexp_by_exposure, glm_estimates, sandwich_estimates )
  lm_out = c(exposure, outcome, n, Wstat, Breusch_Pagan_P, 
             lm_estimates, sandwich_estimates,
             model_R2, model_adjR2, eta_sq, covar_P ) 

  # return(glm_out)
  return(lm_out)

}
