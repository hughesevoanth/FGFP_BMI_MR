################################################################
## A Function to build the AB, trunc_AB, and PA 
##  MT or microbial trait data.
## by: David Hughes
## date: March 16th 2022
################################################################
build_mts = function( wdata , 
                      PA_min_obs = 200,   ## The minimum number of PRESENT of ABSCENT Observations to build a PA trait
                      trunc_AB_present_count = 200,
                      trunc_AB_absence_count = 200,
                      trunc_AB_mean = 20,
                      AB_zero_present_proportion = 0.5,
                      AB_mean_of_zero_trunc_data = 20
                          ){
  
  ## define taxa (column) names
  mt_names = colnames(wdata)
  
  new_mt_data = c()
  for(mt in mt_names){
    trait_counts = wdata[, mt]
    
    ################################
    ## SUMMARY STATISTICS
    ################################
    ### ZERO COUNTS (ABSENCE COUNTS)
    zero_count = sum( trait_counts == 0, na.rm = TRUE )
    
    ## NON-ZERO COUNTS (PRESENCE COUNTS)
    present_count = sum( trait_counts > 0, na.rm = TRUE )
    
    ## PRESENT PROPORTION
    present_proportion = present_count / length(na.omit(trait_counts))
    pp = round(present_proportion, d = 3)
    ## MEAN ABUNDANCE
    trait_mean = mean(trait_counts, na.rm = TRUE)
    
    ## NON-ZERO MEAN ABUNDANCE
    w = which(trait_counts > 0)
    non_zero_trait_mean = mean(trait_counts[w], na.rm = TRUE)
    
    ################################
    ## TRAIT BUILDING
    ################################
    ## PA TRAIT
    if(zero_count >= PA_min_obs & present_count >= PA_min_obs){
      PA = trait_counts
      PA[PA > 1] = 1
      ## name
      n = paste0(mt,"_",pp ,"_PA")
      ## add trait to data set
      new_mt_data = cbind(new_mt_data, PA)
      colnames(new_mt_data)[ncol(new_mt_data)] = n
    }
    
    ## TRUNC AB
    if(present_count >= trunc_AB_present_count & zero_count >= trunc_AB_absence_count & non_zero_trait_mean >= trunc_AB_mean   ){
      trunc_AB = trait_counts
      w = which(  trunc_AB == 0 )
      trunc_AB[w] = NA
      ## name
      n = paste0(mt,"_",pp ,"_trunc_AB")
      ## add trait to data set
      new_mt_data = cbind(new_mt_data, trunc_AB)
      colnames(new_mt_data)[ncol(new_mt_data)] = n
    }
    
    ## AB TRAIT
    if( present_proportion >= AB_zero_present_proportion & non_zero_trait_mean >= AB_mean_of_zero_trunc_data ){
      AB = trait_counts
      # name
      n = paste0(mt,"_",pp ,"_AB")
      ## add trait to data set
      new_mt_data = cbind(new_mt_data, AB)
      colnames(new_mt_data)[ncol(new_mt_data)] = n
    }
  }
  
  return(new_mt_data)
}

