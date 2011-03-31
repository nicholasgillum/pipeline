library(combinat)

setClass("impulse_fit",
         representation(coefs="numeric", mdlerr="numeric", permerr="numeric"))

sigfit <- function(time, expr, dirname) {
  n <- length(time)
  perms <- t(array(unlist(permn(n)),dim=c(n,gamma(n+1))))
  fit <- fitme(time,expr)
  fiterr <-sum( (f(fit,time)-expr)^2)

  d<-dim(perms)
  mdlerr <- matrix(fiterr)
  dir.create(dirname)
  
  for(i in 1:(d[1]/10)) {
    base = (i-1)*10;
    tmp = perms[base:(base+10),]
    start <- Sys.time();
    out<-apply(tmp, 1, fitme_err,time=time);
    write.table(out,paste(dirname,"/","err", i, sep=""))
    end <- Sys.time();
    cat("Total time for round", i," of ",
          d[1]/10, ": ", end-start, "\n")
    }
    mdlerr<-rbind(mdlerr, out);
#  return(new("impulse_fit", coefs=fit, mdlerr=fiterr, permerr=errs));
}

Sig<-function(beta, x, t) {
  return(1/(1+exp(-beta*(x-t))));
}

s1f<-function(p,x) {
  h0=p[1]; h1=p[2]; h2=p[3];
  bt=p[4]; t1=p[5]; t2=p[6];
  return(h0+(h1-h0)*Sig(bt, x, t1));
}

s2f <-function(p,x) {
  h0=p[1]; h1=p[2]; h2=p[3];
  bt=p[4]; t1=p[5]; t2=p[6];
  return(h2+(h1-h2)*Sig(-bt,x,t2))
}

f <- function(p, x){
  h0=p[1]; h1=p[2]; h2=p[3];
  bt=p[4]; t1=p[5]; t2=p[6];
  s1 = s1f(p,x)
  return(1/h1 * s1f(p,x) * s2f(p,x));
}

dfdh0 <- function(p,x) {
  h0=p[1]; h1=p[2]; h2=p[3];
  beta=p[4]; t1=p[5]; t2=p[6]
  s2 = s2f(p,x);
  return( -1/h1 * (1-Sig(beta,x,t1))*s2 );
}

dfdh1 <- function(p,x) {
  h0=p[1]; h1=p[2]; h2=p[3];
  bt=p[4]; t1=p[5]; t2=p[6];
  s1 = s1f(p,x);
  s2 = s2f(p,x);
  return(-(1/h1)*(1/h1)*s1*s2+1/h1*(Sig(bt,x,t1)*s2+s1*Sig(-bt,x,t2)));
}

dfdh2 <- function(p,x) {
  h0=p[1]; h1=p[2]; h2=p[3];
  bt=p[4]; t1=p[5]; t2=p[6];
  s1 = s1f(p,x);
  return(-1/h1 * (1-Sig(-bt,x,t2))*s1);
}

dfdt1 <- function(p,x) {
  h0=p[1]; h1=p[2]; h2=p[3];
  bt=p[4]; t1=p[5]; t2=p[6];
  s2 = s2f(p,x);
  return(1/h1*(-bt*(h1-h0)*Sig(bt,x,t1)*(1-Sig(bt,x,t1)))*s2);
}

dfdt2 <- function(p,x) {
  h0=p[1]; h1=p[2]; h2=p[3];
  bt=p[4]; t1=p[5]; t2=p[6];
  s1 = s1f(p,x);
  return(1/h1*(-bt*(h1-h2)*Sig(-bt,x,t2)*(1-Sig(-bt,x,t2)))*s1);
}

dfdbt <- function(p,x) {
  h0=p[1]; h1=p[2]; h2=p[3];
  bt=p[4]; t1=p[5]; t2=p[6];
  s1 = s1f(p,x);
  s2 = s2f(p,x);
  part1 <- s2/h1*(h1-h0)*(t1-x)*Sig(bt,x,t1)*(1-Sig(bt,x,t1));
  part2 <- s1/h1*(h1-h2)*(t2-x)*Sig(-bt,x,t2)*(1-Sig(-bt,x,t2));
  return(part1+part2);
}

gradf <- function(p,x,signal) {
  stop("Not sure this function works right...")
  differ <- matrix(f(p,x) - signal);
  grad <- rbind(dfdh0(p,x),dfdh1(p,x), dfdh2(p,x),
               dfdt1(p,x),dfdt2(p,x), dfdbt(p,x));
  return(grad %*% differ);
}

objfunc <- function(p, x, signal) {
  mdl <- f(p,x)
  return( 1/2 * sum((mdl - signal)^2) )
}

fitme_err <- function(time, expr) {
  p<-fitme(time,expr)
  return( sum( (f(p,time)-expr)^2 ) )
}

fitme <- function(time, expr) {
  h0 <- expr[1]
  h1 <- max(expr)
  h2 <- expr[length(expr)]
  t1 <- 1/3*(max(time)-min(time))
  t2 <- 2/3*(max(time)-min(time))
  init <-c(h0, h1, h2, rnorm(1), t1, t2)

  control <- list(trace=0,maxit=10^6)
  best<-optim(par=init, fn=objfunc, x=time, signal=expr,
               control=control)

  t=seq(min(time),max(time),length.out=3*length(time))
  interp = approx(x=time, y=expr, xout=t)
  time = interp$x
  expr = interp$y
  
  for(i in 1:100) {
    hg <- runif(3, min=1.5*min(expr), max=1.5*max(expr))
    bg <- rnorm(1, mean=-2)
    tg <-runif(2, min(time), max(time))
    if(tg[1] > tg[2]) { tg<-rev(tg) }
    init <-unlist(c(hg,bg,tg))
    guess <- optim(init,fn=objfunc, x=time, signal=expr,
                   control=control);
    
    if(guess$value < best$value){
      best <- guess;
    }
    
  }
  return(unlist(best$par));
  
}
