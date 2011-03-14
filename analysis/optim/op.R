rm(list=ls(all=TRUE))
library(stats)

mdl_fit<-function(p,x){
  return(p[1]+(p[2]+p[3]*x)*exp(-p[4]*x))
}

func<-function(p,x,y) {
  mdl<-p[1]+(p[2]+p[3]*x)*exp(-p[4]*x)
  return( sum((mdl-y)^2) )
}

x=matrix(unlist(read.table("../time.txt")), ncol=1) #seq(0,5,length.out=npts)
x=x[3:length(x)]
#x-min(x)
y=read.table('../cluster2.txt')
y=matrix(unlist(y[1,3:length(y[1,])]),ncol=1)


fn<-function(p) {
  return(func(p,x,y))
}

p<-c(5, -10.957, -.545,.045)
fit<- optim(p, fn, control=c(trace=1, maxit=10^6), method="SANN")
plot(x,y)
mx<-seq(min(x), max(x), length.out=100);
lines(mx, mdl_fit(fit$par, mx),col='red')

