fhpred_discrete_folds_kai=function(traindata_f,traindata_h,trainlabels,testdata_f,testdata_h,testlabels,preselect_r2,cf,subject)
{
library('MASS')
## down-select to discriminable channels only for feature space-decided not to explicitly include face vs house
rsa = function(d1,d2) {
  # Here is the return arguments list return(outrsa) 
  # function outrsa=rsa(d1,d2) 
  # this function calculates the signed r-squared cross-correlation for 
  # vectors d1 and d2. it is signed to reflect d1>d2 
  
  outrsa <- ((mean(d1)-mean(d2))^3)/abs(mean(d1)-mean(d2))/var(c(d1,d2))*(length(d1)*length(d2))/ (length(d1)+length(d2))^2;
  
  return(outrsa)
} 

ncolr=size(traindata_f[[1]],2)

rf0=matrix(NA, nrow=1,ncol=1)
rh0=matrix(NA, nrow=1,ncol=1)

for  (k in (1:size(traindata_f[[1]],2))){
rf0[k]=rsa(traindata_f[[1]][which(trainlabels==2),k],traindata_f[[1]][which(trainlabels==0),k]);# 2 is face
rh0[k]=rsa(traindata_h[[1]][which(trainlabels==1),k],traindata_h[[1]][which(trainlabels==0),k]);# 1 is house
}

f2u=which(abs(rf0)>preselect_r2);
h2u=which(abs(rh0)>preselect_r2);


## can change this later to be customized classifier if desired
testdata=abind(testdata_f[,f2u],testdata_h[,h2u],along=2)
trainingdata=abind(traindata_f[[1]][which(trainlabels>0),f2u],traindata_h[[1]][which(trainlabels>0),h2u],along=2)
training_label= trainlabels[which(trainlabels>0)];
testing_label= testlabels[which(testlabels>0)];

#t=matrix(rnorm(200*26)/100000,200,26)  
# #training=trainingdata+t
# library('e1071')
# model<-svm(as.factor(training_label)~.,trainingdata)
# res<-predict(model,newdata=testdata,type='decision')
# accuracy=sum(res==testing_label)/length(testing_label)

      
z <- lda(training_label ~ ., as.data.frame(trainingdata))
res<-predict(z,newdata=as.data.frame(testdata),type='decision')
accuracy=sum(res$class==testing_label)/length(testing_label)

# model <- glm(training_label ~.,family=binomial(link='logit'),data=as.data.frame(trainingdata))
#  res<-predict(model,newdata=as.data.frame(testdata),type='response')
#  res=res+1
# accuracy=sum(res==testing_label)/length(testing_label)
out=cbind(res$class, testing_label)
write.csv(out, file = paste("output/",subject,"_label_",cf,".csv",sep=""))
return (accuracy <- accuracy)
       }


