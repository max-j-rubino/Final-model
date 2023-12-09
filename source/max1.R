#max1 fun, makes it so prob is not greater than 1
max1<-function(a,b){
  c=sum(a,b)
  d=min(c,1)
  return(d)
}

