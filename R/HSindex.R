#'
#'The Hopkins and skellam's aggregation index
#'
#'@param sp.ppp an ppp object of a population
#'@param nfocals the number of focal points/events
#'@param k the order of neighbor
#'@param q the number of equal angle sector
#'@param bordercorrection the distance from the border without focal points
#'
#'


HSindex=function(sp.ppp,nfocals,k=1,q=1,bordercorrection=0){
  win=sp.ppp$window
  obs=data.frame(x=sp.ppp$x,y=sp.ppp$y)
  border_xmin=win$xrange[1]+bordercorrection
  border_xmax=win$xrange[2]-bordercorrection
  border_ymin=win$yrange[1]+bordercorrection
  border_ymax=win$yrange[2]-bordercorrection
  
  
  validi=which(sp.ppp$x>border_xmin & sp.ppp$x < border_xmax & 
                 sp.ppp$y>border_ymin & sp.ppp$y< border_ymax)
  if(length(validi)<nfocals){
    etoe_flag=FALSE
    print("number of samples larger than the observed individuals")
    return(NA)
  }
  nqfocals=round(nfocals/q)
  pfocals=data.frame(x=runif(nqfocals,border_xmin,border_xmax),
                     y=runif(nqfocals,border_ymin,border_ymax))
  pr_samples=knq(pfocals,obs,k=k,q=q,type="ptoe")
  
  efocali=sample(validi,nqfocals)
  efocals=data.frame(x=sp.ppp$x[efocali],y=sp.ppp$y[efocali])
  er_samples=knq(efocals,obs,k=k,q=q,type="etoe")
  
  hsindex=sum(pr_samples^2)/sum(er_samples^2)
  return(hsindex)
}
