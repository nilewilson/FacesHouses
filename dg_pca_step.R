dg_pca_step = function(spectra){
  # Here is the return arguments list return(list(pc_weights=pc_weights,pc_vecs=pc_vecs,pc_vals=pc_vals,f=f)) 
  # function [pc_weights,pc_vecs,pc_vals,f]=dg_pca_step(patient,spectra) 
  #this function calculates and returns the principal spectra/eigenvalues,and their projections. 
  nc <- size(spectra,2);#number of channels left 
  ##################################### 
  #create indices to exclude around harmonics of 60 
  f <- 1:300
  no60 <-matrix(data=NA,nrow=0,ncol=1);
  for (k in 1:ceiling(max(f/60)))
  {
    no60 <- append(no60, c((60*k-3):(60*k+3)))
  }
  
  # max(A) takes the max of the whole matrix, If you need column-wise mean, use apply(A,2,max) 
  no60 <- append(no60,c(247:253))
  f <- f[-no60] #dispose of 60hz stuff 
  f=f[-which(f>200)]
  ##################################### 
  ncomps <- length(f);#number of components to keep 
  
  
  #initialize 
  #pc_weights<-zeros(ncomps,nc,size(spectra,3));#projection weights 
  
  pc_vals <- matrix(0, nrow = ncomps, ncol = nc);#eigenvalues 
  pc_vecs<-array(0, dim=c(ncomps,nc,length(f)));#eigenvectors
  
  
  ##run pca 
  
  for (chan in 1:nc){
    tmp1=spectra[f,chan,];
   ts=tmp1[,colSums(abs(tmp1)>0)>1,drop=FALSE]
 #ts=tmp1   
    E=eigen(ts %*% t(ts),TRUE)
    vecs=t(E$vectors)
    vals=t(E$values)
    pc_vecs[,chan,]=vecs[,1:ncomps]
    pc_vals[,chan]=vals[1:ncomps]
  }
  
  return(pc_vecs <- pc_vecs)
} 

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
