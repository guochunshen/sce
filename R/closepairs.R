

closepairs2 <- function(X,rmax,ordered=TRUE){
 
  xy=as.vector(matrix(c(X$x,X$y),nrow=2,ncol=X$n,byrow=TRUE))
  re=nearest_neighbors(xy,rmax)
  names(re)=c("i","j")
  #in c code, i and j is start from 0, thus plus 1 here for consistent index prefer in R
  re$i=re$i+1
  re$j=re$j+1
  re$xi=X$x[re$i]
  re$yi=X$y[re$i]
  re$xj=X$x[re$j]
  re$yj=X$y[re$j]
  re$dx=re$xj-re$xi
  re$dy=re$yj=re$yi
  re$d=sqrt((re$xj-re$xi)^2+(re$yj-re$yi)^2)
  if(!ordered){
    ok=re$i<re$j
    re$i=re$i[ok]
    re$j=re$j[ok]
    re$xi=re$xi[ok]
    re$yi=re$yi[ok]
    re$xj=re$xj[ok]
    re$yj=re$yj[ok]
    re$dx=re$dx[ok]
    re$dy=re$dy[ok]
    re$d=re$d[ok]
  }
  return(re)
}