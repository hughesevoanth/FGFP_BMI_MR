################################################################
## A Function to derive a few summary statistics and filter
##  16S count data for those to retain for analysis.
## by: David Hughes
## date: March 16th 2022
################################################################
taxa_2_keep = function(wdata , ## data frame with MT count data
                       taxa_col_names , ## the column names to interegate
                       non_zero_count = 100 , ## the number of samples with non-zero count data
                       zero_trunc_mean = 25 , ## the mean abundance across the zero-truncated abundance data
                       single_individual_min = NA ## the minimum count seen in at least 1 sample to retain taxa
){
  
  ## Summary Statistics
  sumstats = t(sapply(taxa_col_names, function(mt){
    x = na.omit( study_data$mt_data[, mt] )
    cv = var(x,na.rm=TRUE)/mean(x, na.rm = TRUE)
    ## Identify the non-zero observations
    w = which(x > 0)
    x = x[w]
    ## Non-Zero count
    count = length(x)
    ## min, mean, max
    min = min(x)
    mean = mean(x)
    max = max(x)
    trunc_cv = var(x,na.rm=TRUE)/mean(x, na.rm = TRUE)
    ##
    out = c(count, min, mean, max, cv, trunc_cv)
    names(out) = c("count", "min", "mean", "max", "cv", "trunc_cv")
    return(out)
  }) )
  ## Make sum stats a data frame
  sumstats = data.frame(sumstats)
  
  ## Taxa to keep
  if(!is.na(single_individual_min)){
    w = which( sumstats$count >= non_zero_count & (sumstats$mean >= zero_trunc_mean | sumstats$max > single_individual_min) )  
  } else {
    w = which( sumstats$count >= non_zero_count & sumstats$mean >= zero_trunc_mean  )  
  }
  
  ## Taxa (column) names to keep
  out = taxa_col_names[w]
  
  return(out)
  
}

