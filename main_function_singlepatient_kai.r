main_function_singlepatient_kai=function(subject){

library(R.matlab)
dataset<-readMat(paste("data/", subject, "_faceshouses.mat", sep=""))
  
stim=dataset$stim;
srate=dataset$srate;
data=dataset$data
source('car.R')
data=car(data)
data=array(unlist(data), dim = c(nrow(data[[1]]), ncol(data[[1]])))
################raw signal##########################
# source('main_function_ws_modified_kai.R')
# #source('main_function_ws_pre.R') is for Kai's code
# output=main_function_ws_modified_kai(data,stim)

##########
source('main_function_ws_modified.R')
output=main_function_ws_modified(data,stim)

source('fhpred_discrete_folds_kai.R')
# accuracy_raw=0;
# for (cf in 1:3)
# {
#   accuracy_raw[cf]=fhpred_discrete_folds_kai(traindata_f=output$f_template_train_fold[cf],
#                                              traindata_h=output$h_template_train_fold[cf],
#                                              
#                                              trainlabels=output$train_events_fold[[cf]][,2],
#                                              
#                                              testdata_f=output$f_template_test_fold[[cf]],
#                                              testdata_h=output$h_template_test_fold[[cf]],
#                                              testlabels=output$test_events_fold[[cf]][,2],
#                                              preselect_r2=.05)
# }
# 
# avg_accuracy_raw=mean(accuracy_raw)



# #############bb signal####################
source('fh_get_events.R')
source('calc_dg_spectra_updated.R')
source('calc_nspectra.R')
source('dg_pca_step.R')
source('kjm_lnA_timecourse.R')
source ('fh_pc_clean.R')
# 
pts=fh_get_events(stim)
pts=array(unlist(pts), dim = c(nrow(pts[[1]]), ncol(pts[[1]])))
if (pts[1,1]<500){ pts=pts[-1,]}

spectra=calc_dg_spectra_updated(data,pts) ##
spectra=array(unlist(spectra), dim = c(nrow(spectra[[1]]), ncol(spectra[[1]]),size(spectra[[1]],3)))

##test
#  sp_m<-readMat('spectra.mat')
#  spectra_m=sp_m$spectra

nspectra=calc_nspectra(spectra);##ok
nspectra=array(unlist(nspectra), dim = c(nrow(nspectra[[1]]), ncol(nspectra[[1]]),size(nspectra[[1]],3)))
 
# ns<-readMat('nspectra.mat')
# nspectra_m=ns$nspectra
#  
 
 pc_vecs=dg_pca_step(nspectra) ##ok
#  ns<-readMat('pc_vecs.mat')
# pc_vecs_m=ns$pc.vecs
#  
# ## get broadband timecourse from 1st spectral principle component
# 
# # create indices to exclude around harmonics of 60Hzplot
f0 <- 1:200
no60 <-matrix(data=NA,nrow=0,ncol=1);
for (k in 1:ceiling(max(f0/60)))
{
  no60 <- append(no60, c((60*k-3):(60*k+3)))
}

# max(A) takes the max of the whole matrix, If you need column-wise mean, use apply(A,2,max) 
no60 <- append(no60,c(247:253))
f0 <- f0[-no60] #dispose of 60hz stuff 

# get lnA
lnA=matrix(0,nrow=size(data,1),ncol=size(data,2))
for (chan in 1:size(data,2))
{
dt=data[,chan];
tmp1=pc_vecs[,chan,]
mm=t(tmp1[,colSums(abs(tmp1)>0)>1,drop=FALSE])
#mm=tmp1
pcvec1=mm[,1]
##for test ##

lnA[,chan]=kjm_lnA_timecourse(dt,pcvec1,srate,f0);# this generates lnA timeseries for each channel
}
# 
# # smooth lnA,exponentiate,and subtract 1 to get bb timecourse
# lnA2<-readMat('lnA.mat')
# lnA_m=lnA2$lnA
# 

bb=fh_pc_clean(-lnA);
 bb=bb-1;


###call preprocessing function
# b<-readMat('bb.mat')
# bb_m=b$bb

output_bb=main_function_ws_modified(bb,stim)

####make prediction###
##load raw data features###
f_template_test_fold_raw=output$f_template_test_fold
h_template_test_fold_raw=output$h_template_test_fold

f_template_train_fold_raw=output$f_template_train_fold
h_template_train_fold_raw=output$h_template_train_fold

train_events_fold=output$train_events_fold
test_events_fold=output$test_events_fold

###load bb data features
f_template_test_fold_bb=output_bb$f_template_test_fold
h_template_test_fold_bb=output_bb$h_template_test_fold

f_template_train_fold_bb=output_bb$f_template_train_fold
h_template_train_fold_bb=output_bb$h_template_train_fold

##

##create lists to save final results
f_template_test_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))
h_template_test_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))
f_template_train_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))
h_template_train_fold<-lapply(1:3, function(x) matrix(NA, nrow=2, ncol=3))

for (cf in 1:3){
f_template_test_fold[[cf]]=abind(f_template_test_fold_raw[[cf]],f_template_test_fold_bb[[cf]],along=2)
h_template_test_fold[[cf]]=abind(h_template_test_fold_raw[[cf]],h_template_test_fold_bb[[cf]],along=2)

f_template_train_fold[[cf]]=abind(f_template_train_fold_raw[[cf]],f_template_train_fold_bb[[cf]],along=2)
h_template_train_fold[[cf]]=abind(h_template_train_fold_raw[[cf]],h_template_train_fold_bb[[cf]],along=2)

}
  
  accuracy_bth=0;
for (cf in 1:3)
{
  accuracy_bth[cf]=fhpred_discrete_folds_kai(traindata_f=f_template_train_fold[cf],
                                     traindata_h=h_template_train_fold[cf],
                                     
                                     trainlabels=train_events_fold[[cf]][,2],
                                     
                                     testdata_f=f_template_test_fold[[cf]],
                                     testdata_h=h_template_test_fold[[cf]],
                                     testlabels=test_events_fold[[cf]][,2],
                                     preselect_r2=.05, cf=cf, subject=subject)
}

#avg_accuracy_bth=mean(accuracy_bth)
#avg_accuracy=append(avg_accuracy_raw,avg_accuracy_bth)
# 
 write.csv(accuracy_bth, file = paste("output/", subject,"_accuracy_bth.csv",sep=""))
 file=paste("output/output_",subject,".RData")
 save.image(file)
rm(list=ls()) 
}
