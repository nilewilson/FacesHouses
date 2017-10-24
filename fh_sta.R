fh_sta = function(inputdata,events,fh_class,tlims) {
  
  library(PopED)
  cls_times= events[which(events[,2]==fh_class),1]
  
  sta=matrix(data=0, nrow=(tlims[2]-tlims[1]+1),ncol=size(data,2))
  
  for (k in 1:length(cls_times)){
  sta=sta+inputdata[cls_times[k]+c(tlims[1]:tlims[2]),];
}
  
  sta=sta/k;
  
return(sta)
} 
