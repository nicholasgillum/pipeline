library(stats)
source('mdl.R')
data<-read.table("exp_data.txt")
time<-data$time-min(data$time)
expr<-data[,2]#data$VNG1365C#



fn<-function(p,x,y) {
  return(sum( (log_mdl(p,x) - y)^2 ))
}

start<-optim(c(max(expr), min(expr)*(min(expr)>0),
               .0,0,2), fn, x=time,y=expr,
             control=c(trace=1,maxit=10^6),method="SANN")
fit<-optim(start$par, fn, x=time,y=expr,
           control=c(trace=1,maxit=10^4))

mtime<-seq(min(time),max(time), length.out=1000)
plot(time,expr,col="black")
lines(mtime, log_mdl(fit$par, mtime))
