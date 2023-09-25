sys <- Sys.info()[1]

if (sys=="Windows") {
 setwd("..\\Slicker\\examples\\R_example")  #Set path to input file directory 
}

target <- read.delim("test5_input_target.txt",header=FALSE,sep=" ")      #2D array
timebase <- read.delim("test5_input_timebase.txt",header=FALSE)          #1D reconstruction times
p1 <- read.delim("test5_input_proxy1.txt",header=FALSE,sep=" ")          #2D proxy 
p2 <- read.delim("test5_input_proxy2.txt",header=FALSE,sep=" ")          #2D proxy 

source("../SLICKER.R", keep.source=TRUE)   # Load and run SLICKER function

recon <- SLICKER(target,timebase,p1,p2,num_ensemble=512,time_limit=20,tol=1e-7)   #Exceute Win or UNIX reconstruct

write.table(recon,"test5_R_output.txt",sep="\t",row.names=FALSE,quote=FALSE,
col.names=c("!time base","M-estimator","+/- for_95%_CI","Qn"))
