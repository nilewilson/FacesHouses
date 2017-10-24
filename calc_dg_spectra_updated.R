

calc_dg_spectra_updated = function(data,pts)  {
  # Here is the return arguments list return(list(spectra=spectra)) 
  #this function calculates the the spectra at time points in pts(:,2) 
  
  ########################################################################### 
  #parameters 
  library(R.matlab)
  library(e1071)
  library('bspec')
  #library(multitaper)
  samplefreq <- 1000;#sampling frequency 
  bsize <- 1000;#window size for spectral calculation 
  ########################################################################## 
  #calculate spectra 
  #print('calculating power spectra',quote = FALSE) 
  spectra<-array(0, dim=c(300,size(data,2),size(pts,1))) 
  
  for (k in 1:size(pts,1)){ 
    #if (mod[k,100]==0),sprintf(1,'%03d ',k);if (mod[k,500]==0),sprintf(1,'*/%d\r',size[pts,1]);end,end 
    for (m in 1:size(data,2)){
     tmp=as.ts(data[(pts[k,2]-floor(bsize/2)+1):(pts[k,2] + ceiling(bsize/2)),m]);
     ts<- welchPSD(tmp, seglength=800) ## is different than psd in matlab
     power=ts$power*(ts$frequency)/2
     spectra[,m,k] <- as.vector(power[2:301]) ## why 2:301?
     
#       tt=data[(pts[k,2]-floor(bsize/2)+1):(pts[k,2] + ceiling(bsize/2)),m];
#       pgram01 <- spec.pgram(ts(tt, frequency=samplefreq), plot=FALSE)
#       
#     
#       spectra[,m,k] <- as.vector(pgram01$spec[2:301])
    } 
  }
  
  ########################################################################## 
  return(list(spectra <- spectra))
} 
