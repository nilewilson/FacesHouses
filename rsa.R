rsa = function(d1,d2) {
# Here is the return arguments list return(outrsa) 
# function outrsa=rsa(d1,d2) 
# this function calculates the signed r-squared cross-correlation for 
# vectors d1 and d2. it is signed to reflect d1>d2 

outrsa <- ((mean(d1)-mean(d2))^3)/abs(mean(d1)-mean(d2))/var(c(d1,d2))*(length(d1)*length(d2))/ (length(d1)+length(d2))^2;

return(outrsa)
} 

