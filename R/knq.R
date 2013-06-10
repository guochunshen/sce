#'
#' To find the kth nearest neighborhood in the q equal angle sector for each forcal point/event
#' 
#'@param focals a data.frame with colnames x and y that contains spaital locaiton of each focals
#'@param obs the observed population distribution, a data.frame with colnames x and y
#'@param k the order of neighborhood
#'@param q the number of equal angle sectors
#'@param type the type of neighborhood distance. it can be "ptoe" (point to event distance) or "etoe" (event to event distance)
#'
#'



knq=function(focals,obs,k=1,q=1,type="ptoe"){

  n=dim(obs)[1]
  m=dim(focals)[1]
  focalx=rep(focals$x,each=n)
  focaly=rep(focals$y,each=n)
  d=sqrt((focalx-obs$x)^2+(focaly-obs$y)^2)
  selfi=d==0
  d[selfi]=NA
  
  dim(d)=c(n,m)
  re=matrix(nrow=m,ncol=q)
  if(q==1){
    re[,1]=apply(d,2,function(x) (sort(x))[k])
  }else{
    angle=atan2(obs$y-focaly,obs$x-focalx)
    angle[angle<0]=angle[angle<0]+2*pi
    breaks=seq(0,2*pi,length.out=q+1)
    dbks=as.numeric(cut(angle,breaks,include.lowest=TRUE))
    #dbks[selfi]=NA
    dim(dbks)=c(n,m)
    
    for(i in 1:m){
      for(j in 1:q){
        re[i,j]=sort((d[,i])[dbks[,i]==j])[k] 
      }
    }
  }
  return(re)
}
