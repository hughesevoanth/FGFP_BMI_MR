################################################################
## A Function to rank normal transform a vector
##  
## by: David Hughes
## date: March 16th 2022
################################################################
rntransform = function(x , split_ties = TRUE){
  if (split_ties == TRUE) {
    out <- rank(x, ties.method = "random") - 0.5
  }
  else {
    out <- rank(x) - 0.5
  }
  out[is.na(x)] <- NA
  out <- out/(max(out, na.rm = T) + 0.5)
  out <- qnorm(out)
  out
}

