
#' 
#' calcualte a pair wise distance matrix and stored in a filebacked big matrix
#'
#'
#'
dist_bigmatrix<-function(x,y,filename){
  n=length(x)
  A=big.matrix(n,n,type="double",backingfile=filename)
  for(i in 1:n){
    A[i,]=sqrt((x[i]-x)^2+(y[i]-y)^2)
  }
  return(A)
}
