#'
#'Fitting the each group individual mapped population in the community by a clustering model
#'
#'\code{fitCluster} fit the point pattern of each groups (species, DBH class) by a homogeneous/heterogeneous Cox/cluster process.
#' 
#'
#'@param com a \link{scp} object of a given community
#'@param trend a specific formular for modeling trend.
#'@param cluster define how internal cluster assembled with intensity function. possible choose are "poisson" and "LGCP"
#'@param sigTest a logical flag to test significance of extra aggregation in residual and habitat variables
#'@param ctlpars control parameters used in model fitting. See meaning of each par in Details.
#'
#'
#'@details
#'Three steps were carried out:
#'1. fit the point pattern by a heterogeneous poisson process; 2. calculate pair correlation function of the 
#' residual. 3. estimate internal clustering parameters by minimum contrast method.
#'
#'@return a list of fm object contains fitted parameters for each group
#'
#'@seealso 
#'\link{sigTest}, \link{removeRareSpecies}
#'
#'
#'@examples
#' #load the testData set
#' data(testData)
#' #remove rare species
#' com=removeRareSpecies(testData,80)
#' 
#' #fit pattern of each species by a best cluster model
#' fittedmodel=fitCluster(com,~elev+grad)
#' 
#'

fitCluster<-function(com,trend=~1,cluster="LGCP",sigTest=FALSE,
                     ctlpars=list("rmax"=25,"rmin"=3,"bw"=5,"sigma2"=3,"alpha"=10,
                                  "nurange"=c(Inf,0.5),"q"=2,"edgecor"='translate',
                                  "nsim"=10,"r"=seq(0,60,2)),...){
  #validation check
  if(!RandomFieldsSafe()){
    stop("The newest version of RandomFields package is needed")
  }
  #number of environmental variables
  nhabitat=length(com$habitat)
  
  data.ppp=com$com
  #fit the heterogeneous poisson point process
  data.ppm=ppm(data.ppp,trend,covariates=com$habitat)
  
  beta=coef(data.ppm)
  
  if("poisson"!=cluster){
    #estimate the intensity of points at each point location
    lambda=predict(data.ppm,locations=data.ppp, type="trend")
    #estimate pair correlation function by adaptive method
    g=pcf_adaptive(data.ppp,maxt=ctlpars$rmax,lambda,bw=ctlpars$bw,adaptive=0.5,kerneltype=1)
    estPars=best.matern.estpcf(g, c(sigma2=ctlpars$sigma2, alpha=ctlpars$alpha),
                               rmax=ctlpars$rmax,rmin=ctlpars$rmin,nu=ctlpars$nurange,
                               q=ctlpars$q)
    nu=attr(estPars,"nu")
    names(nu)="nu"
    sigma2=estPars$par[1]
    alpha=estPars$par[2]
    re=c(nu,sigma2,alpha,beta)
    attr(re,"minicontrast")=estPars
  }else{
    re=c(rep(NA,3),beta)
  }
  attr(re,"data")=com
  attr(re,"trend")=trend
  attr(re,"fittedmodel")=data.ppm
  attr(re,"modeltype")=cluster
  attr(re,"ctlpars")=ctlpars
  attr(re,"class")=c("fm",class(re))
  if(sigTest){
    aggreRes_pvalue=sigAggreResidualTest(re,ctlpars$nsim,ctlpars$r,ctlpars$edgecor)
    names(aggreRes_pvalue)="aggreRes"
    clusterResidual= (aggreRes_pvalue<ctlpars$siglevel) & cluster!="poisson"
    habitat_pvalues=sigHabitatTest(re,clusterResidual)
    attr(re,"pvalues")=c(aggreRes_pvalue,habitat_pvalues)
  }
  
  return(re)
}

best.matern.estpcf=function(X, startpar = c(sigma2 = 1, alpha = 1), lambda = NULL, nu=c(0.25), 
                            q = 1/4, p = 2, rmin = NULL, rmax = NULL, ...){
  bestmodel=list()
  minivalues=numeric()
  for(i in 1:length(nu)){
    bestmodel[[i]]=matern.estpcf(X, startpar,rmax=rmax,rmin=rmin,nu=nu[i],q=q)
    minivalues[i]=bestmodel[[i]]$opt$value
  }
  minione=which(minivalues==min(minivalues))[1]
  return(bestmodel[[minione]])
}

matern.estpcf=function (emppcf, startpar = c(sigma2 = 1, alpha = 1), lambda = NULL, nu=1/2, 
                        q = 1/4, p = 2, rmin = NULL, rmax = NULL, ...) 
{
  
  startpar <- check.named.vector(startpar, c("sigma2", "alpha"))
  
  theoretpcf <- function(par,r,nu){
    if (any(par <= 0) | par[2]>5000) 
      return(rep(Inf, length(r)))
    
    if(nu==Inf){
      return(exp(par[1]*exp(-1/2*(r/par[2])^2)))
    }else{
      return(exp(Covariance(r,model="matern",param=c(0.0,par[1],0.0,par[2],nu))))
    }
  }
  
  result <-  mincontrast(emppcf, theoretpcf, startpar, 
                         ctrl = list(q = q, p = p, rmin = rmin, rmax = rmax), 
                         fvlab = list(label = "%s[fit](r)", desc = "minimum contrast fit of LGCP"), 
                         ...,nu=nu)
  
  par <- result$par
  if(FALSE){
    plot(emppcf)
    theo=theoretpcf(par,emppcf$r,nu)
    lines(x=emppcf$r,y=theo,col=3)
  }
  names(par) <- c("sigma2", "alpha")
  result$par <- par
  mu <- if (is.numeric(lambda) && length(lambda) == 1 && lambda > 
              0) 
    log(lambda) - par[["sigma2"]]/2
  else NA
  result$modelpar <- c(sigma2 = par[["sigma2"]], alpha = par[["alpha"]], 
                       mu = mu)
  result$internal <- list(model = "matern")
  attr(result,"nu")=nu
  return(result)
}

print.fm<-function(x,...){
  xname=names(x)
  x=as.numeric(x)
  names(x)=xname
  print(x)
}

