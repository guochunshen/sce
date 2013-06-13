#' To calculate pvalues of habitat varaibles under consideration of extra aggregation existed in residual.
#' 
#' @param fittedModel a \link{fm} object representing a fitted cluster point process for the given point pattern.
#' @param clusterResidual a logical flag to consider any left cluster residual if it is true.
#'
#'
#'@details
#' if there is significant clustered residual (\code{clusterResidual=TRUE}), significance of the habitat will be
#' evaludated by considering those extra aggregation in residual. otherwise (\code{clusterResidual=FALSE}), regular
#' significant test method of regression coefficients was carried out.
#'
#'
#'
#'@examples
#'
#'data(testData)
#'
#'sp1=subset(testData,testData$traits$species=="ACALDI")
#'
#'fm=fitCluster(sp1,~elev+grad,sigTest=FALSE)
#'
#'sigHabitatTest(fm)
#'

sigHabitatTest<-function(fittedModel,clusterResidual=TRUE){
  if(as.character(attr(fittedModel,"trend"))[2]=="1"){
    warning("can't test significant of habitat without included it in the model")
    pvalue=1
    names(pvalue)="habitat"
    return(pvalue)
  }
    data.ppm=attr(fittedModel,"fittedmodel")
    nu=fittedModel[1]
    alpha=fittedModel[3]
    sigma2=fittedModel[2]
    beta=fittedModel[-c(1:3)]
    if(clusterResidual)
      acacov=vcov.mykppm(data.ppm,par=c(sigma2,alpha),nu=nu)
    else
      acacov=vcov(data.ppm)
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
  # evaluate integral
  ff <- Z * lambda * wt
  J <- t(Z) %*% ff
  
  # compute pair correlation function minus 1
  E=pcfoneless (U, par, nu, nU, ff) 
  
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

pcfoneless <- function (U, par, nu, nU, ff,filename="temp_backingfile",usingBigmatrix=FALSE) {
  #since the bigalgebra still not avaiable on the linux,just throw an error when handle very large nU
  #remove it once these package are avalibled in Linux again
  canBigmatrix=require(bigmemory) & require(bigalgebra)
  if(canBigmatrix){
    nUlimit=5000
  }else{
    nUlimit=800000
  }
  
  if((!usingBigmatrix) & nU<nUlimit){
    r <- as.vector(pairdist(U))
    gr <- g(r,par,nu)-1
    rm(r)
    gc()
    G <- matrix(gr, nU, nU)
    rm(gr)
    gc()
    E <- t(ff) %*% G %*% ff
    rm(G)
    gc()
  }else{
    if(!canBigmatrix){
      stop("too large number of menmory need in test significance of habitat")
    }
    
    neach=1e3
    n=floor(nU/neach)
    #in case of parallal calculation, backfile name should not be the same for different cores,
    #Thus add some random number to the filename
    filename=paste(filename,as.character(runif(1)),sep="")
    des_filename=paste(filename,".desc",sep="")
    G=big.matrix(nU,nU,type="double",backingfile=filename,descriptorfile=des_filename)
    x=rep(U$x,each=neach)
    y=rep(U$y,each=neach)
    for(i in 1:n){
      sel=((i-1)*neach+1):(neach*i)
      G[sel,]=g(sqrt((U$x[sel]-x)^2+(U$y[sel]-y)^2),par,nu)-1
    }
    if(nU>n*neach){
      sel=(n*neach+1):nU
      x=rep(U$x,each=length(sel))
      y=rep(U$y,each=length(sel))
      G[sel,]=g(sqrt((U$x[sel]-x)^2+(U$y[sel]-y)^2),par,nu)-1
    }
    E <- (t(ff) %*% G %*% ff)[,]
    #removing the filebacking object and file
    rm(G)
    gc()
    file.remove(filename)
    file.remove(des_filename)
  }
  
  return(E)
}
# 
# 
# gloop=function(r,par,nu){
#   nr=length(r)
#   nmax=1e8
#   if(nr>nmax){
#     neach=nmax
#     n=floor(nr/neach)
#     re=numeric()
#     for(i in 1:n){
#       sel=((i-1)*neach+1):(neach*i)
#       re=c(re,g(r[sel],par,nu))
#     }
#     if(nr>n*neach){
#       sel=(n*neach+1):nr
#       re=c(re,g(r[sel],par,nu))
#     }
#     return(re)
#   }else{
#     return(g(r,par,nu))
#   }
#   
# }

g <- function(r, par,nu) {
  
  if(nu==Inf){
    return(exp(par[1]*exp(-1/2*(r/par[2])^2)))
  }else
    return(exp(Covariance(r,model="matern",param=c(0.0,par[1],0.0,par[2],nu))))
  
}

