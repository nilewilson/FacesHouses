kjm_lnA_timecourse=function(data,pcvec1,srate,freqs){
  # function tf=kjm_tf_pwr(data,srate,freqs)generate complex time-frequency
  # analysis in range freqs
  #
  
  ##
  library(fftw)
  ##for test ##
#   freqs=f0;
#  data=dt;
#   
  tfp=zeros(length(data),length(freqs));#time-frequency power
  data_1=data;
  for (tmp in 1:length(freqs))
   # for (tmp in 1:20)
  { 
    freq=freqs[tmp];
    ## create morlet wavelet
    t=1:floor(5*srate/freq);
    tmid=floor(max(t)/2);
    
    wvlt_1=exp(1i*2*pi*(freq/srate)*(t-tmid))*exp(-(t-tmid)^2/(2*(srate/freq)^2));#gaussian envelope
    
    ## calculate convolution
    wvlt=append(wvlt_1,zeros(1,(length(data_1)-1)))
    
    data=append(data_1, zeros(1,(length(wvlt_1)-1)))
    tconv=IFFT(FFT(wvlt)*FFT(data))
    ############################################
    
    v1=(1:(floor(length(wvlt_1)/2)-1));
    v2=(floor(length(tconv)-length(wvlt_1)/2+1):length(tconv));
    tconv=tconv[-v2]
    tconv=tconv[-v1]
    
    tconv=abs(tconv)^2;# power
    #if (mean(tconv)==0) break #if there is some problem
    tconv=tconv/mean(tconv);#normalize power at this freq 
    
    tfp[,tmp]=tconv;
    
  }
  
  pc1=log(tfp)%*%pcvec1 
  return(pc1<-pc1)
}