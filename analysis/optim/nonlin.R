exp_mdl <- function(time, A, B, R) {
  return(A+B*exp(R*time))
}

gen_exp<-function(time,A=2,B=10,R=-1/2) {
  perfect<- exp_mdl(time,A,B,R)
  noise <- rnorm(length(time), 0, 1/20*(max(perfect)-min(perfect)))
  return(perfect+noise)
}

crit_mdl<-function(time, A,B,C,R) {
  return(A+(B+C*time)*exp(R*time));
}

gen_crit <- function(time, A=2,B=1,C=2, R=-1/2) {
  perfect<-crit_mdl(time, A,B,C,R);
  noise<-rnorm(length(time), 0, 1/20*(max(perfect)-min(perfect)))
  return(perfect+noise)
}

fit_expr<-function(time, expr) {
  return(lm(log(expr)~log(time),x=TRUE, y=TRUE))
}

fit_crit<-function(time, expr) {
  return(lm(log(expr)~(log(time)*time),x=TRUE, y=TRUE))
}

plot_fit <- function(time, expr, lmfit) {
  plot(time,expr)
  points(time, fitted(lmfit),col="red")
}
