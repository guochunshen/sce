#'
#'transfer a marked ppp object into a scp object, it mark is just the species name
#'
#'@param data.ppp an ppp object with species name as its mark
#'

as.scp<-function(data.ppp){
  if(is.null(data.ppp$marks)){
    sp=rep(1,data.ppp$n)
  }else{
    if(is.data.frame(data.ppp$marks))
      sp=as.vector(data.ppp$marks[,1])
    else
      sp=as.vector(data.ppp$marks)
  }
  
  return(scp(species=sp,x=data.ppp$x,y=data.ppp$y,win=data.ppp$window))
}