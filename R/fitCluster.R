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
#'1. fit the point pattern by a heterogeneous poisson process; 
#'2. calculate pair correlation function of the residual (using adaptive estimation, see details in \code{\link{pcf_adaptive}}). 
#'3. estimate internal clustering parameters by minimum contrast method.
#'
#'the \code{ctlpars} contains parameter values that are very sensitive to the final result. it controls calculations of 
#'the pari correlation function, the minimum contrast method and significance tests of habitat and extra internal aggregation.
#'it designed for advanced user, thus change it when you know what you are doing.
#'
#'Specifically, \code{rmax} controls the spatial range that calculate the pair correlation function in the
#' \code{\link{pcf_adaptive}} function. \code{rmax} and \code{rmin} also define the range of distance considered in the 
#' minimum contrast method n estimating parameters of internal aggregation. \code{nurange} controls the sharp of marten 
#' covariance function in the minimum control method. \code{bw} controls the bindwidth of kernel estimation
#' in \code{\link{pcf_adaptive}}. \code{sigma2}, \code{alpha} and \code{edgecor} are initial parameters of 
#' the fitted clustering process.  \code{nsim}, \code{r}, \code{ntry} and \code{edgecor} are parameters used in testing significance of aggregated
#' pattern in the residual (see \code{\link{sigAggreResidualTest}}). it adapted a goodness-of-fit test method. Thus they control the number of simulations, interpoint
#' ditance range and edge correction method respectively. 
#' 
#' aditionally, \code{nd} controls the maximum dimension of the grid of dummy points (nd * nd or nd[1] * nd[2]) used to evaluate the integral
#'  in the pseudolikelihood. Incompatible with eps. without control it, you might have memory limitation problem in the \code{\link{sigHabitatTest}} 
#'  function. 
#'
#'@return a list of fm object contains fitted parameters for each group
#'
#'@seealso 
#'\link{sigHabitatTest}, \link{sigAggreResidualTest},\link{updateCluster} ,\link{varDecomp}
#'
#'
#'@examples
#' #load the testData set
#' data(testData)
#' 
#' #remove rare species
#' com=removeRareSpecies(testData,600)
#' 
#' #fit pattern of each species by a best cluster model
#' fittedmodel=fitCluster(com,~elev+grad,sigTest=TRUE)
#' 
#'

fitCluster<-function(com,trend=~1,cluster="LGCP",sigTest=FALSE,
                     ctlpars=list("rmax"=25,"rmin"=1,"bw"=2,"sigma2"=3,"alpha"=10,
                                  "nurange"=c(Inf,0.5),"edgecor"='translate',
                                  "nsim"=10,"r"=seq(0,60,1),"siglevel"=0.05,
                                  nd=50,ntry=10),...){
  #validation check
  if(!RandomFieldsSafe()){
    stop("The newest version of RandomFields package is needed")
  }
  #number of environmental variables
  nhabitat=length(com$habitat)
  
  data.ppp=com$com
  
  #try to avoide memory problem in sigHabitatTest function.
  dQ=default.dummy(data.ppp)
  if(dQ$n>10000){
    #fit the heterogeneous poisson point process
    data.ppm=ppm(data.ppp,trend,covariates=com$habitat,nd=ctlpars$nd)
  }else{
    #fit the heterogeneous poisson point process
    data.ppm=ppm(data.ppp,trend,covariates=com$habitat)
  }
  
  
  beta=coef(data.ppm)
  
  if("poisson"!=cluster){
    #estimate the intensity of points at each point location
    lambda=predict(data.ppm,locations=data.ppp, type="trend")
    #estimate pair correlation function by adaptive method
    g=pcf_adaptive(data.ppp,maxt=ctlpars$rmax,lambda,bw=ctlpars$bw,adaptive=0.5,kerneltype=1)
    estPars=best.matern.estpcf(g, c(sigma2=ctlpars$sigma2, alpha=ctlpars$alpha),
                               rmax=ctlpars$rmax,rmin=ctlpars$rmin,nu=ctlpars$nurange,
                               q=2)
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
    #some very strange distribution will cause some error in the simulation
    attr(re,"pvalues")=sigTestofCluster(re)
  }
  
  return(re)
}

sigTestofCluster <- function (re) {
  
  ctlpars=attr(re,"ctlpars")
  cluster=attr(re,"modeltype")
  aggreRes_pvalue=try(sigAggreResidualTest(re,ctlpars$nsim,ctlpars$r,ctlpars$edgecor,ntry=ctlpars$ntry))
  clusterResidual= (aggreRes_pvalue<ctlpars$siglevel) & cluster!="poisson"
  habitat_pvalues=try(sigHabitatTest(re,clusterResidual))
  if(inherits(aggreRes_pvalue,"try-error") | inherits(habitat_pvalues,"try-error")){
    attr(clusterpvalues,"possible reason")="extreme unusual spatial distribution"
    warning("There are errors in calculation of pvalues, it might be caused by extreme unusual spatial distribution")
    re=c(NA,NA)
    class(re)="try-error"
    return(re)
  }else{
    return(c(aggreRes_pvalue,habitat_pvalues))
  }
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

