main_function_ws_modified=function(data,stim){
  
  library("R.matlab")
#  source('fhpred_discrete_folds.R')
  library(abind)

##create lists to save final results
f_template_test_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))
h_template_test_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))
f_template_train_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))
h_template_train_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))
train_events_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))
test_events_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))
sta_f_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))
sta_h_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))
stim_train_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))
stim_test_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))

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
  #evs<= evs[which(evs[ ,1]<500|evs[ ,2]>(length(stim)-1000)), ] ##not surese
  return(evs)
} 

## create events vector 
events=fh_get_events(stim);
events= events[which(events[,3]!=0),]
events=events[,-2]

signal=data;
srate=1000;
tlength=nrow(signal);
num_chans=ncol(signal);
tlims=c(-199,400);# times to start and end erps
erp_baseline=c(-199,50);# times to calcualte erp based upon(must be within tlims)

## identify thirds

# first third
ind1=1:floor(tlength/3);
ev1=(events[,1]<ind1[length(ind1)])&(events[,1]>ind1[1]);
ev1=events[ev1,];

# second third
ind2=(floor(tlength/3)+1):floor(tlength*2/3);
ev2=((events[,1]<ind2[length(ind2)])&(events[,1]>ind2[1]));
ev2=events[ev2,];

#third third
ind3=(floor(tlength*2/3)+1):tlength;
ev3=((events[,1]<ind3[length(ind3)])&(events[,1]>ind3[1]));
ev3=events[ev3,]

## 3 times cross-folding-divide into train and test segments



for (cf in 1:3) {
  if (cf==1)
  {
    # get test event indices
    test_t=ind3;
    train_t=c(ind1,ind2);
    test_e=ev3;
    train_e=abind(ev1,ev2,along=1);
  }
  else if (cf==2) 
  {# get test event indices
    test_t=ind2;
    train_t=append(ind1,ind3);
    test_e=ev2;
    train_e=abind(ev1,ev3,along=1);
  }
  else if (cf==3) 
  {# get test event indices
    test_t=ind1;
    train_t=append(ind2,ind3);
    test_e=ev1;
    train_e=abind(ev2,ev3,along=1)
  } 
  
  
  
  # scale by std,etc(based only on train segments)
  for (k in 1:num_chans){
    signal[,k]=signal[,k]/sd(signal[train_t,k]);# signal(:,k)=(signal[:,k]-mean[signal[train_t,k]])/std(signal[train_t,k]);
  }
  
  # get class specific{face,house}STA templates from train-zero out for 100 ms pre-stim
  ##function of fh_sta
  
  fh_sta = function(inputdata,events,fh_class,tlims) {
    cls_times= events[which(events[,2]==fh_class),1]
    sta=matrix(data=0, nrow=(tlims[2]-tlims[1]+1),ncol=ncol(inputdata))
    for (k in 1:length(cls_times)){
      sta=sta+inputdata[cls_times[k]+c(tlims[1]:tlims[2]),];
    }
    sta=sta/k;
    return(sta)
  } 
  
  
  
  #get sta templates
  ##tlims=[-199 400] for debug, it should be the result from the main func
  tlims=c(-199, 400)
  sta_h=fh_sta(signal,train_e,1,tlims);# houses
  sta_f=fh_sta(signal,train_e,2,tlims);# faces
  
  # recenter stas w.r.t. baseline
  for (k in 1:num_chans) {
    sta_h[,k]=sta_h[,k]-mean(sta_h[c(erp_baseline[1]:erp_baseline[2])-tlims[1]+1,k]);
    sta_f[,k]=sta_f[,k]-mean(sta_f[c(erp_baseline[1]:erp_baseline[2])-tlims[1]+1,k]);
  }
  
  ## generate train data-
  # get samples of template at appropriate point and at every 100ms for each class-
  # note that this method is not efficient,because templates are class specific. 
  
  # generate 4 training points between each stimulus,
  # must be min 100ms from any actual stimulus point and min 50 from each other
  a=train_e[,1];
  a=sort(a);
  npts=matrix(data=NA,nrow=0,ncol=1);
  library(PopED)
  library(pracma)
  for (k in 2:length(a)){
    b=randperm(floor((a[k]-a[k-1]-200)/100));# floor(rand[4,1]*[a[k]-a[k-1]-200])
    b=100*b[1:4]+floor(50*rand(4,1))+a[k-1];
    npts=rbind(npts,b);
  }
  
  train_e=abind(train_e, abind(npts,0*npts,along=2),along=1);
  
  # get projection into templates
  f_template_train=zeros(size(train_e,1),num_chans);
  h_template_train=0*f_template_train;
  #
  for (k in 1:size(train_e,1)){
    # dot products
    for (chan in 1:num_chans){
      dt=signal[train_e[k,1]+c(tlims[1]:tlims[2]),chan];# select data
      dt=dt-mean(dt[c(erp_baseline[1]:erp_baseline[2])-tlims[1]+1]);#baseline data
      
      f_template_train[k,chan]=sum(sta_f[,chan]*dt);# convolve
      h_template_train[k,chan]=sum(sta_h[,chan]*dt);# convolve
    }}
  
  ## generate test data-note that loop method sucks in matlab
  
  f_template_test=zeros(size(test_e,1),num_chans);
  h_template_test=0*f_template_test;
  
  
  for (k in 1:size(f_template_test,1)){# note that these times line up with a phase lag according to tlims(1)
    for (chan in 1:num_chans)
    {
   # dt=signal[k+tlims[1]:tlims[2],chan];#select data (wrong)
    dt=signal[test_e[k,1]+c(tlims[1]:tlims[2]),chan]
    dt=dt-mean(dt[c(erp_baseline[1]:erp_baseline[2])-tlims[1]+1]);#baseline data
    f_template_test[k,chan]=sum(sta_f[,chan]*dt);#convolve
    h_template_test[k,chan]=sum(sta_h[,chan]*dt);# convolve
    }}
  
  ## store in appropriate fold,
  
  
  f_template_test_fold[cf]=list(f_template_test);
  h_template_test_fold[cf]=list(h_template_test);
  #
  f_template_train_fold[cf]=list(f_template_train);
  h_template_train_fold[cf]=list(h_template_train);
  #
  train_events_fold[cf]=list(train_e);
  test_events_fold[cf]=list(test_e);
  #
  sta_f_fold[cf]=list(sta_f);
  sta_h_fold[cf]=list(sta_h);
  #
  stim_train_fold[cf]=list(stim[train_t]);
  stim_test_fold[cf]=list(stim[test_t[1:size(f_template_test,1)]]);
}
newlist<-list ("f_template_test_fold"=f_template_test_fold,"h_template_test_fold"=h_template_test_fold,
               "f_template_train_fold"=f_template_train_fold,"h_template_train_fold"=h_template_train_fold,
               "train_events_fold"=train_events_fold,"test_events_fold"=test_events_fold,
               "sta_f_fold"=sta_f_fold,"sta_h_fold"=sta_h_fold,"stim_train_fold"=stim_train_fold,
               "stim_test_fold"=stim_test_fold)
return(newlist<-newlist)

}

## save testing and training data for later classification 
#save('data/' subject '/' subject '_' cls '_cross_folds','*fold*','\n') 
#write.csv(f_template_test_fold,'f_template_test.csv');
#h_template_test_fold
#
#f_template_train_fold[cf]=list(f_template_train);
#h_template_train_fold[cf]=list(h_template_train);
#
#train_events_fold[cf]=list(train_e);
#test_events_fold[cf]=list(test_e);
