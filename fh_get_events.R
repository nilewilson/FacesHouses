fh_get_events = function(stim)  {
# Here is the return arguments list return(list(evs=evs)) 
#this function defines evs lying at the beginning(1st column)and middle 
#(2nd column)of each stimulus and isi period 
# 0=ISI 
# 1=HOUSE STIM 
# 2=FACE STIM 
  ##
tmp=c(0,stim[1:((length(stim)-1))])
b <- which((stim-tmp)!=0 ) 
c <- floor(diff(b)/2) 
b<-b[1:length(b)-1]  
d <- b + c 
evs=matrix(data=NA, nrow=length(b),ncol=3)

evs[ ,1] <- b 
evs[ ,2] <- d 
evs[ ,3] <- stim[d] 
evs<-evs[which( evs[ ,3] != 0 ), ]
evs[which( evs[ ,3] < 51 ),3] <- 1 
evs[which( evs[ ,3] == 101 ),3] <- 0 
evs[which( evs[ ,3] > 50 ),3] <- 2 
rm(b,c,d) 
 
#clip if too close to ends of run files 
# evs<= evs[which(evs[ ,1]<500|evs[ ,2]>(length(stim)-1000)), ] ##not sure
return(list(evs <- evs))
} 
