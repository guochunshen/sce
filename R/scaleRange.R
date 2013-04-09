#'
#'scale a vector's range and keep the relationship between elements
#'
#'@param x a numeric vector or maxtrix to be scale the range
#'@param min the minimum value of the output vector
#'@param max the maxmum value of the output vector
#'
#'
#'@examples
#'
#'x=runif(100,1,100)
#'x2=scaleRange(x,0.5,1)
#'cor(x,x2)
#'
#'x=matrix(x,10,10)
#'x3=scaleRange(x,0.5,1)
#'cor(as.numeric(x),as.numeric(x3))



scaleRange<-function(x,min,max){
  gmax=max(x)
  gmin=min(x)
  a=(gmax-gmin)/(max-min)
  b=(gmin*max-gmax*min)/(max-min)
  exp_x=(x-b)/a
  return(exp_x)
}
