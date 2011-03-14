library(splines)
library(stats)
library(boot)

brokenstick<-function(time, expr) {
  lmfit <- lm(expr ~ time)

  below<-2:(length(time)-3)
  above<-length(time)-below-1
  texpr<-expr[below+1,]
  ttime<-time[below+1,]
  rsqrd<-matrix(rep(0,length(ttime)),ncol=1)
  for( i in 1:length(ttime) ) {
    pivot<-ttime[i]
    fm1<-lm(expr ~ time, subset=(time<=pivot),x=TRUE,y=TRUE)
    fm2<-lm(expr ~ time, subset=(time>pivot), x=TRUE,y=TRUE)
    rs1<-corr(cbind(matrix(fm1$y,ncol=1), matrix(fitted(fm1),ncol=1)))
    rs2<-corr(cbind(matrix(fm2$y,ncol=1), matrix(fitted(fm2),ncol=1)))
    rsqrd[i,]=(above[i]+below[i])^-1 * (rs1*(below[i]+1) + rs2*above[i])
  }
  
  return(ttime[rsqrd==max(rsqrd)])
}
par(mfrow=c(2,1))
expr_data<-read.table('cluster2.txt')
time<- matrix(unlist(read.table('time.txt')),ncol=1);
expr<-matrix(unlist(expr_data[1,]),ncol=1)
pivot <- brokenstick(time, expr)
fm1 <- lm(expr~time, subset=(time <= pivot))
fm2 <- lm(expr~time, subset=(time >  pivot))
plot(time,expr)
abline(coefficients(fm1))
abline(coefficients(fm2))

time <- matrix(seq(1,100,length.out=100),ncol=1)
expr <- matrix(5*time*(time<=50) + (2*time+150)*(time>50)+rnorm(100,0,5),ncol=1)
pivot<-brokenstick(time,expr)

plot(time,expr)
fm1 <- lm(expr~time, subset=(time <= pivot))
fm2 <- lm(expr~time, subset=(time >  pivot))

abline(coefficients(fm1))
abline(coefficients(fm2))

