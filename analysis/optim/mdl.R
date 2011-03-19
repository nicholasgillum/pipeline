glr_mdl <-function(p,x) {
	return(p[1]/(1+exp(p[2]+p[3]*x)));
}

exp_mdl<-function(p,x) {
  return(p[1]*(1-exp(p[2]*(x-p[3])))*(x>=p[3]))
}

log_mdl<-function(p,x) {
  return(p[1]+(p[2]-p[1])/( (1-exp(-p[3]*(x-p[4])))^(1/p[5]) ))
}
