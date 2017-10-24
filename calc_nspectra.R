calc_nspectra = function(spectra){
  # Here is the return arguments list return(list(spectra=spectra)) 
  #this function normalizes the spectra prior to the pca step 
  #kjm 12/07 
  #spectra=readMat('spectra.mat')
  library(R.matlab)
  library(pracma)
  # normalization step 
  for (k in 1:size(spectra,2)){
    
    single_spectra=spectra[,k,]
    m=rowMeans(single_spectra,2);
   
    
    #tmp=m[,colSums(m>0)>1,drop=FALSE]
    
#    t2=repmat(m,1,size(single_spectra,2))
    
    
    t2=matrix(repmat(m,1,size(single_spectra,2)), nrow=size(single_spectra,1),ncol=size(single_spectra,2))
    
    spectra[,k,]=log(single_spectra/t2);
  } 
  rm(k) 
  return(list(spectra <- spectra))
} 
###############################
