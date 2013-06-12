#'Generate simulated population with negative binomial abundance distribution
#'
#'@param lam expected individual intensity of a population
#'@param a the aggregated parameter in the negative binomial distribution. the smaller the a is, the more aggregate the populaiton is.
#'@param win a owin object represent the range of the population
#'@param nxy a vector of number of quadrats in the x and y axis respectively.
#'
#'
#'@return
#'a scp object
#'
#'

rPopNeg=function(lam,a,win,nxy){
  n=nxy[1]*nxy[2]
  xbreaks=seq(win$xrange[1],win$xrange[2],length.out=nxy[1]+1)
  ybreaks=seq(win$yrange[1],win$yrange[2],length.out=nxy[2]+1)
  gridarea=(xbreaks[2]-xbreaks[1])*(ybreaks[2]-ybreaks[1])
  
  neach=rnbinom(n,mu=lam*gridarea,size=1/a)
  x=numeric()
  y=numeric()
  for(i in 1:nxy[1]){
    for(j in 1:nxy[2]){
      ij=j+(i-1)*nxy[1]
      x=c(x,runif(neach[ij],xbreaks[i],xbreaks[i+1]))
      y=c(y,runif(neach[ij],ybreaks[j],ybreaks[j+1]))
    }
  }
  re=scp(species=rep(1,length(x)),x=x,y=y,win=win)
  return(re)
}