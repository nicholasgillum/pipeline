library(stats)
library(boot)
source('mdl.R')

dataframe <- read.table("exp_data.txt")

critfit<-function(dataframe, i) {
  rawtime <- dataframe$time/10
  rawdata <- dataframe[,i]

  premask <- rawtime <=0;

  ptime <- rawtime[premask];
  pexpr <- rawdata[premask];

  time  <- rawtime[!premask]
  expr  <- rawdata[!premask]
}

glrfit <- function(dataframe, i) {
  rawtime <- dataframe$time/10
  rawdata <- dataframe[,i]

  premask <- rawtime <= 0;

  ptime <- rawtime[premask];
  pexpr <- rawdata[premask];

  time  <- rawtime[!premask];
  expr  <- rawdata[!premask];
  
  b1 <- max(expr);
  b2 <- 2
  b3 = -1
  
  fn <- function(p,x,y) {
    return(sum( (glr_mdl(p,x)-y) ^ 2 ))
  }
  
  fit<-optim(c(b1,b2,b3), fn, x=time, y=expr,
             control=c(trace=0, maxit=10^4))

  
  
  output = list(ptime=ptime, pexpr=pexpr, time=time,
    expr=expr, glrfit=fit, glrRsq=corr(cbind(expr, glr_mdl(fit$par,time))))
  plot(time, expr, col="red")
  points(time, glr_mdl(fit$par, time),pch=22)
  return(output)
}



