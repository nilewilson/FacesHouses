fh_pc_clean = function(lnA)  {
  # Here is the return arguments list return(bb) 
  #this function smooths the pcs, 
  
  #lnA=readMat(('lnA.mat'))
  library(PopED)
  library(signal)
  ## ##DO SMOOTHING####### 
  winlength <- 80 
  
  bb <-matrix(0, nrow=size(lnA,1),ncol=size(lnA,2))
  for (k in 1:size(lnA,2))
  {
    # convolve with gaussian window 
    q=lnA[,k]
    lnAS=conv(q,gausswin(80))
    # clip edges 
    lnAS=lnAS[-(1:(floor(winlength/2-1)))] 
    lnAS=lnAS[-((length(lnAS)-floor(winlength/2-1)):length(lnAS))]
    
    #z-score 
    lnAs <- (lnAS-mean(lnAS))/sd(lnAS)
    # For robustness c() added to divide, Assume you divide by a scalar, not matrix
    # mean(A) takes the average of the whole matrix, If you need column-wise mean, use colMeans(A) 
    # re-exponentiate(as originally in log) 
    bb[ ,k] <- exp(lnAs) 
  } 
  return(bb)
} 
