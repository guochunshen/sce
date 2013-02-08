#' To calculate pvalues of habitat varaibles under consideration of extra aggregation existed in residual.
#' 
#' @param fittedModel a \link{fm} object representing a fitted cluster point process for the given point pattern.
#'
#'

sigHabitatTest<-function(fittedModel){
  
  data.ppm=attr(fittedModel,"fittedmodel")
  nu=fittedModel[1]
  alpha=fittedModel[2]
  sigma2=fittedModel[3]
  beta=fittedModel[-c(1:3)]
  acacov=vcov.mykppm(data.ppm,par=c(sigma2,alpha),nu=nu)
  pvalue=2*(1-pnorm(abs(beta/sqrt(diag(acacov)))))[-1]
  return(pvalue)
}

vcov.mykppm <- function(object, ...,
                      what=c("vcov"),par,nu)
{
  what <- match.arg(what)
  #verifyclass(object, "kppm")
  # extract composite likelihood results
  po <- object
  # extract quadrature scheme information
  Q <- quad.ppm(po)
  U <- union.quad(Q)
  nU <- npoints(U)
  wt <- w.quad(Q)
  # compute fitted intensity values
  lambda <- fitted(po, type="lambda")
  # extract covariate values
  Z <- model.matrix(po)
  # compute pair correlation function minus 1
  r <- as.vector(pairdist(U))
  gr <- gloop(r,par,nu)-1
  G <- matrix(gr, nU, nU)
  # evaluate integral
  ff <- Z * lambda * wt
  J <- t(Z) %*% ff
  E <- t(ff) %*% G %*% ff
  # asymptotic covariance matrix in the Poisson case
  J.inv <- try(solve(J))
  # could be singular
  if(inherits(J.inv, "try-error")) {
    if(what == "internals") {
      return(list(ff=ff, J=J, E=E, J.inv=NULL))
    } else {
      return(NULL)
    }
  }
  # asymptotic covariance matrix in the clustered case
  vc <- J.inv + J.inv %*% E %*% J.inv
  
  switch(what,
         vcov={ return(vc) },
         corr={
           sd <- sqrt(diag(vc))
           co <- vc/outer(sd, sd, "*")
           return(co)
         },
         fisher={
           fish <- try(solve(vc))
           if(inherits(fish, "try-error")) fish <- NULL
           return(fish)
         },
         internals={
           return(list(ff=ff, J=J, E=E, J.inv=J.inv, vc=vc))
         })
  stop(paste("Unrecognised option: what=", what))
}

gloop=function(r,par,nu){
  nr=length(r)
  nmax=1e8
  if(nr>nmax){
    neach=nmax
    n=floor(nr/neach)
    re=numeric()
    for(i in 1:n){
      sel=((i-1)*neach+1):(neach*i)
      re=c(re,g(r[sel],par,nu))
    }
    if(nr>n*neach){
      sel=(n*neach+1):nr
      re=c(re,g(r[sel],par,nu))
    }
    return(re)
  }else{
    return(g(r,par,nu))
  }
  
}

g <- function(r, par,nu) {
  
  if(nu==Inf){
    return(exp(par[1]*exp(-1/2*(r/par[2])^2)))
  }else
    return(exp(Covariance(r,model="matern",param=c(0.0,par[1],0.0,par[2],nu))))
  
}

