car = function(data)  {
# Here is the return arguments list return(list(data=data)) 
#this function calculates and returns the common avg reference of 2-d matrix "data". 
#kjm 12/07 
 
#data <- double(data) 
 
transflag <- 0 
if (size(data,1) < size(data,2)) { 
    data <- t(data);
    transflag <- 1 ;} 
 
num_chans <- size(data,2)
 
# create a CAR spatial filter 
spatfiltmatrix <- -matrix(1, nrow = num_chans, ncol = num_chans) 

for (i in (1:num_chans)) { 
    spatfiltmatrix[i,i] <- num_chans-1 
} 

spatfiltmatrix <- spatfiltmatrix/num_chans 
 
# perform spatial filtering 
if (all(is.na(spatfiltmatrix))==FALSE )  { 
 #  sprintf(1,'Spatial filtering ') 
   data <- data %*% spatfiltmatrix 
  
    if ( size(data,2) != size(spatfiltmatrix,1))  { 
      sprintf(1,'The first dimension in the spatial filter matrix has to equal the second dimension in the data') 
   } 
} 
 
if (transflag == 1)  { 
data <- t(data) 
} 
return(list(data <- data))
} 
