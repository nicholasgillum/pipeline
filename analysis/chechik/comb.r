generate_orderings<-function(n){
  output <- matrix(rep(n,n), ncol=1);
  for(i in rev(1:(n-1))) {
    output<-rbind(output, matrix(rep(i, i),ncol=1))
  }
  return(output);
}

a<-generate_orderings(3)
