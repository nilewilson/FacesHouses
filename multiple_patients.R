
library(R.matlab)
library(PopED)
library(fftw)

# Set your working directory here
setwd('E:/FacesHouses')

source('main_function_singlepatient_kai.R')

# Run the main function per subject
main_function_singlepatient_kai('de')
main_function_singlepatient_kai('fp')
main_function_singlepatient_kai('ja')
main_function_singlepatient_kai('wc')
main_function_singlepatient_kai('zt')
main_function_singlepatient_kai('mv')
main_function_singlepatient_kai('aa')
main_function_singlepatient_kai('ap')
main_function_singlepatient_kai('rn')
main_function_singlepatient_kai('ha')

