rm(list=ls())
library(stats)
source("defs.r")

gen_data<-function(n) {
  time_data  <- c(-40,-20,0,10,20,40,60,80)/10
  #time_data <- c(0,15,30,45,60,120)
  time_mask  <- time_data >= 0;
  time       <- time_data[time_mask];
  expr_data  <- read.table('ura3_A.txt')
  #expr_data <- read.table('environ.txt')

  starttime <- Sys.time()
  expr <- as.vector(expr_data[n,time_mask],mode="numeric")
  sigfit(time, expr, paste(n))
  endtime <- Sys.time()
  print(endtime-starttime)
}

