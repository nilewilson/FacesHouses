
library(R.matlab)
library(PopED)
library(fftw)

# Set your working directory here
setwd('.')

source('main_function_singlepatient_kai.R')

# Run the main function per subject
main_function_singlepatient_kai('de')
print("de complete")
main_function_singlepatient_kai('fp')
print("fp complete")
main_function_singlepatient_kai('ja')
print("ja complete")
main_function_singlepatient_kai('wc')
print("wc complete")
main_function_singlepatient_kai('zt')
print("zt complete")
main_function_singlepatient_kai('mv')
print("mv complete")
main_function_singlepatient_kai('aa')
print("aa complete")
main_function_singlepatient_kai('ap')
print("ap complete")
main_function_singlepatient_kai('rn')
print("rn complete")
main_function_singlepatient_kai('ha')
print("ha complete")

