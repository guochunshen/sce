# #'fit the given point pattern by a best cluster model 
# #'and then estimate the proportion of variance explained by habitat and internal cluster
# #'
# #' \code{varDecomp} using three steps to estimate proportion of variance in spatial point pattern. First,
# #' it fit the given point pattern by a cluster model, in which significance of each habitat variables,
# #' signficance of remaining cluster were estimated. Then a best cluster model was chosen by the backward
# #' stepwise selection method.
# #'
# #'@param com a \link{scp} object that contains distribution of trees and environmental variables
# #'@param multicore a logical flag indicate whether using multicore to speedup calculation.
# #'        physical support of multicore computation is needed first
# #'@param mc.cores number of cores used. it also means how many species were under calculation simutatnously.
# #'@param pars sensitive internal parameters related to model fitting and significance test.
# #'        do not change it unless what you are doing.
# #'
# #'
# #'
# #'
# 
# 
# 
# 
# #
# # Main functions:
# #  
# # backwardStep(): model selection.
# # varDecomp.per(): estimate parameters for heterogeneous Poisson or Heterogeneous Cox point process.
# # myagregativeResidualTest(): significant test of habitat.
# # vcov.kppm(): calculate p-values for each regression parameters.
# # 
# #
# # Author: guochun shen
# # Email: shenguochun@gmail.com
# ###############################################################################
# 
# 
# varDecomp<-function(com,multicore=FALSE,mc.cores=1,
#                     pars=list("rmax"=70,"rmin"=5,"bw"=15,"nsim"=39,"siglevel"=0.05,
#                                   "sigma2"=3,"alpha"=10,"nurange"=c(Inf,0.5),"q"=2,
#                                   "siglevel2"=0.05,"agtest"=TRUE,"edgecor"='translate')){
#   dolist=1:com$S
#   dolist=sample(dolist)
#   covr=??
#   
#   if(multicore){
#     result=mclapply(dolist,onesp,data=com,covr=covr,sp=com$sp,win=com$win,mc.cores=mc.cores)
#   }else{
#     result=lapply(dolist,onesp,data=com,covr=covr,sp=com$sp,win=com$win)
#   }
#   result=unlist(result)
#   dim(result)=c(length(result)/length(dolist),length(dolist))
#   return(result)
# }
# 
# onesp=function(x,data,covr,sp,win){
#   print(paste("species",x))
#   subdata=data[data$sp==sp[x],]
#   data.ppp=unique(ppp(x=subdata$gx,y=subdata$gy,window=win))
#   re=backwardStep(data.ppp,covr,pars=list("rmax"=70,
#                                           "rmin"=5,"bw"=15,"nsim"=39,"siglevel"=0.05,
#                                           "sigma2"=3,"alpha"=10,"nurange"=c(Inf,0.5),"q"=2,
#                                           "siglevel2"=0.05,"agtest"=TRUE,"edgecor"='translate'))
#   result=c(re,attr(re,"beta"))
#   result=c("pForH"=attr(re,"pForH"),"pForD"=attr(re,"pForD"),result)
#   write.table(result,paste("./temp/tp",x,".txt",sep=""))
#   return(result)
# }
# 
# 
# 
# backwardStep=function(data.ppp,covrs,pars=list("rmax"=60,"rmin"=5,"bw"=10,
#                                                "sigma2"=3,"alpha"=10,"nurange"=c(Inf,0.5),"q"=0.5,
#                                                "nsim"=39,"siglevel"=0.01,"siglevel2"=0.05,
#                                                "agtest"=FALSE,"edgecor"='translate'),confi=FALSE){
#   start=TRUE
#   del=NA
#   cname=names(covrs)
#   plotdim=c(covrs[[1]]$xrange[2],covrs[[1]]$yrange[2])
#   dQ=default.dummy(data.ppp)
#   #avoide memory problem
#   if(dQ$n>10000){
#     dn=100
#   }else{
#     dn=NA
#   }
#   while(start | !is.na(del)){
#     start=FALSE
#     fittedModel=varDecomp.per(data.ppp=data.ppp,covrs=covrs,pars=pars,dn=dn)
#     pvalues=attr(fittedModel,"pvalues")
#     insig=pvalues>pars$siglevel
#     
#     if(any(insig)){
#       del=which(pvalues==max(pvalues,na.rm=T))[1]
#       covrs=covrs[-del]
#       
#     }else{
#       del=NA
#     }
#   }
#   
#   #calculate global pvalue for covariable
#   if(length(covrs)!=0 ){
#     pForH=1-any(pvalues<(1-(1-0.05)^(1/length(cname))))
#     beta=attr(fittedModel,"beta")[match(cname,names(covrs))]
#     beta[is.na(beta)]=0
#     attr(fittedModel,"beta")=beta
#     
#     if(confi){
#       conf.inter=conf.inter(fittedModel,covrs,plotdim,nsim=200,pars)
#       attr(fittedModel,"conf")=conf.inter
#     }
#     
#   }else{
#     pForH=1
#     attr(fittedModel,"beta")=rep(0,length(cname))
#   }
#   
#   attr(fittedModel,"pForH")=pForH
#   return(fittedModel)
# }
# 
# 
# #estimate best parameters for a given data.ppp and covrs
# varDecomp.per=function(data.ppp,covrs,pars,covrsig=TRUE,dn=NA){
#   if(length(covrs)==0){
#     data.ppm=ppm(data.ppp)
#     trend=NA
#     tildeZ=0
#   }else{
#     formu=formu(names(covrs))
#     
#     if(!is.na(dn))
#       data.ppm=ppm(data.ppp,nd=dn,formu,covariates=covrs)
#     else
#       data.ppm=ppm(data.ppp,formu,covariates=covrs)
#     trend=predict(data.ppm, type="trend",ngrid=dim(covrs[[1]])/2)
#     if(any(chgs<-trend$v>10)){
#       trend$v[chgs]=10
#     }
#     tildeZ=var(log(as.numeric(trend$v)))
#   }
#   beta=coef(data.ppm)
#   lambda=predict(data.ppm,locations=data.ppp, type="trend")
#   
#   if(pars$agtest){
#     #browser()
#     pvalueForD=myagregativeResidualTest(data.ppp,covrs,trend,
#                                         pars$nsim,edcor=pars$edgecor)
#   }else
#     pvalueForD=0
#   
#   if(pvalueForD<=pars$siglevel2){
#     #lambda=predict(data.ppm,locations=data.ppp, type="trend")
#     g=mygest(data.ppp,maxt=pars$rmax,lambda,bw=pars$bw,adaptive=0.5,kerneltype=1)
#     K=Kinhom(data.ppp,lambda,correction=pars$edgecor,nlarge=5000,r=c(0,g$t))
#     K=K[-1,]
#     g=data.frame(r=g$t,theo=rep(1,length(g$t)),iso=g$g)
#     attributes(g)=attributes(K)
#     #estPars=best.matern.estBoth(K, g, c(sigma2=pars$sigma2, alpha=pars$alpha),
#     #		rmax=pars$rmax,rmin=pars$rmin,nu=pars$nurange,q=1)
#     estPars=best.matern.estpcf(g, c(sigma2=pars$sigma2, alpha=pars$alpha),
#                                rmax=pars$rmax,rmin=pars$rmin,nu=pars$nurange,q=pars$q)
#     nu=attr(estPars,"nu")
#     sigma2=estPars$par[1]
#     if(!covrsig){
#       return(c(tildeZ,sigma2))
#     }
#     #K=Kinhom(data.ppp,lambda,correction=pars$edgecor,nlarge=5000,r=seq(0,120,1))
#     alpha=matern.estK(K,c(alpha=pars$alpha),sigma2=sigma2,rmax=pars$rmax,
#                       rmin=pars$rmin,nu=nu,q=0.25)$par
#     print(c(estPars$par[2],alpha))
#     #browser()
#     acacov=vcov.kppm(data.ppm,par=c(sigma2,alpha),nu=nu)
#     re=c("tildeZ"=tildeZ,sigma2,alpha,"nu"=nu)
#   }else{
#     acacov=vcov(data.ppm)
#     re=c("tildeZ"=tildeZ,NA,NA,NA)
#   }
#   
#   #calculate pvalue
#   pval=2*(1-pnorm(abs(beta/sqrt(diag(acacov)))))[-1]
#   attr(re,"pvalues")=pval
#   attr(re,"pForD")=pvalueForD
#   attr(re,"beta")=beta[-1]
#   return(re)
# }
# 
# 
# 
# myagregativeResidualTest=function(data.ppp,covrs,trend,nsim,edcor){
#   
#   if(length(covrs)!=0){
#     e=envelope(data.ppp,fun=pcf,simulate=expression(rpoispp(trend)),
#                savefuns=TRUE,r=seq(0,60,2),nsim=nsim,verbose=FALSE,correction=edcor)
#   }else{
#     e=envelope(data.ppp,fun=pcf,verbose=FALSE,
#                savefuns=TRUE,r=seq(0,60,2),nsim=nsim,correction=edcor)
#   }
#   npp=nsim+1
#   Kfuns=attr(e,"simfuns")
#   n=length(e$obs)
#   Kfuns[,1]=e$obs
#   ui=vector("numeric",npp)
#   for (k in 2:(n-1)){
#     Kit=Kfuns[k,]
#     t.rslt=step.ui(Kit,2,npp)
#     ui = ui+t.rslt
#   }
#   pvalue=calc.pval(ui)
#   return(pvalue)
# }
# 
# 
# matern.estpcf=function (emppcf, startpar = c(sigma2 = 1, alpha = 1), lambda = NULL, nu=1/2, 
#                         q = 1/4, p = 2, rmin = NULL, rmax = NULL, ...) 
# {
#   
#   startpar <- check.named.vector(startpar, c("sigma2", "alpha"))
#   
#   theoretpcf <- function(par,r,nu){
#     if (any(par <= 0) | par[2]>5000) 
#       return(rep(Inf, length(r)))
#     
#     if(nu==Inf){
#       return(exp(par[1]*exp(-1/2*(r/par[2])^2)))
#     }else{
#       return(exp(Covariance(r,model="matern",param=c(0.0,par[1],0.0,par[2],nu))))
#     }
#   }
#   
#   result <-  mincontrast(emppcf, theoretpcf, startpar, ctrl = list(q = q, 
#                                                                    p = p, rmin = rmin, rmax = rmax), fvlab = list(label = "%s[fit](r)", 
#                                                                                                                   desc = "minimum contrast fit of LGCP"), 
#                          ...,nu=nu)
#   
#   par <- result$par
#   if(FALSE){
#     plot(emppcf)
#     theo=theoretpcf(par,emppcf$r,nu)
#     lines(x=emppcf$r,y=theo,col=3)
#   }
#   names(par) <- c("sigma2", "alpha")
#   result$par <- par
#   mu <- if (is.numeric(lambda) && length(lambda) == 1 && lambda > 
#               0) 
#     log(lambda) - par[["sigma2"]]/2
#   else NA
#   result$modelpar <- c(sigma2 = par[["sigma2"]], alpha = par[["alpha"]], 
#                        mu = mu)
#   result$internal <- list(model = "matern")
#   attr(result,"nu")=nu
#   return(result)
# }
# 
# 
# best.matern.estpcf=function(X, startpar = c(sigma2 = 1, alpha = 1), lambda = NULL, nu=c(0.25), 
#                             q = 1/4, p = 2, rmin = NULL, rmax = NULL, ...){
#   bestmodel=list()
#   minivalues=numeric()
#   for(i in 1:length(nu)){
#     bestmodel[[i]]=matern.estpcf(X, startpar,rmax=rmax,rmin=rmin,nu=nu[i],q=q)
#     minivalues[i]=bestmodel[[i]]$opt$value
#   }
#   minione=which(minivalues==min(minivalues))[1]
#   return(bestmodel[[minione]])
# }
# 
# 
# matern.estK=function (X, startpar = c(alpha = 1),sigma2 = 1, lambda = NULL, nu=1/2, 
#                       q = 1/4, p = 2, rmin = NULL, rmax = NULL, ...) 
# {
#   dataname <- deparse(substitute(X))
#   
#   startpar <- check.named.vector(startpar, c("alpha"))
#   Integrand <- function(r, par,nu,sigma2) {
#     #here is the modification
#     if(nu==Inf){
#       return(2*pi*r*exp(sigma2*exp(-1/2*(r/par[1])^2)))
#     }else{
#       return(2*pi*r*exp(Covariance(r,model="matern",param=c(0.0,sigma2,0.0,par[1],nu))))
#     }
#   }
#   theoret <- function(par, rvals, ..., integrand,nu,sigma2) {
#     if (any(par <= 0) | par[1]>5000) 
#       return(rep(Inf, length(rvals)))
#     th <- numeric(length(rvals))
#     th[1] <- if (rvals[1] == 0) 
#       0
#     else integrate(integrand, lower = 0, upper = rvals[1], 
#                    par = par,nu=nu,sigma2=sigma2)$value
#     for (i in 2:length(rvals)) th[i] = th[i - 1] + integrate(integrand, 
#                                                              lower = rvals[i - 1], upper = rvals[i], par = par,
#                                                              nu=nu,sigma2=sigma2)$value
#     return(th)
#   }
#   result <-  mincontrast(X, theoret, startpar, ctrl = list(q = q, 
#                                                            p = p, rmin = rmin, rmax = rmax), fvlab = list(label = "%s[fit](r)", 
#                                                                                                           desc = "minimum contrast fit of LGCP"), explain = list(dataname = dataname, 
#                                                                                                                                                                  fname = attr(X, "fname"), modelname = "Cox process with matern pair correlation function"), 
#                          ..., integrand = Integrand,nu=nu,sigma2=sigma2)
#   
#   par <- result$par
#   names(par) <- c("alpha")
#   result$par <- par
#   
#   return(result)
# }
# 
# 
# 
# vcov.kppm <- function(object, ...,
#                       what=c("vcov"),par,nu)
# {
#   what <- match.arg(what)
#   #verifyclass(object, "kppm")
#   # extract composite likelihood results
#   po <- object
#   # extract quadrature scheme information
#   Q <- quad.ppm(po)
#   U <- union.quad(Q)
#   nU <- npoints(U)
#   wt <- w.quad(Q)
#   # compute fitted intensity values
#   lambda <- fitted(po, type="lambda")
#   # extract covariate values
#   Z <- model.matrix(po)
#   # compute pair correlation function minus 1
#   r <- as.vector(pairdist(U))
#   gr <- gloop(r,par,nu)-1
#   G <- matrix(gr, nU, nU)
#   # evaluate integral
#   ff <- Z * lambda * wt
#   J <- t(Z) %*% ff
#   E <- t(ff) %*% G %*% ff
#   # asymptotic covariance matrix in the Poisson case
#   J.inv <- try(solve(J))
#   # could be singular
#   if(inherits(J.inv, "try-error")) {
#     if(what == "internals") {
#       return(list(ff=ff, J=J, E=E, J.inv=NULL))
#     } else {
#       return(NULL)
#     }
#   }
#   # asymptotic covariance matrix in the clustered case
#   vc <- J.inv + J.inv %*% E %*% J.inv
#   
#   switch(what,
#          vcov={ return(vc) },
#          corr={
#            sd <- sqrt(diag(vc))
#            co <- vc/outer(sd, sd, "*")
#            return(co)
#          },
#          fisher={
#            fish <- try(solve(vc))
#            if(inherits(fish, "try-error")) fish <- NULL
#            return(fish)
#          },
#          internals={
#            return(list(ff=ff, J=J, E=E, J.inv=J.inv, vc=vc))
#          })
#   stop(paste("Unrecognised option: what=", what))
# }
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
# 
# g <- function(r, par,nu) {
#   
#   if(nu==Inf){
#     return(exp(par[1]*exp(-1/2*(r/par[2])^2)))
#   }else
#     return(exp(Covariance(r,model="matern",param=c(0.0,par[1],0.0,par[2],nu))))
#   
# }
# 
# 
# rLGCP=function(en.filter,sigma2,alpha,N,plotdim,nu=0.5){
#   sigma2true=sigma2
#   covrdim=dim(en.filter$v)
#   xcell=plotdim[1]/covrdim[2]
#   ycell=plotdim[2]/covrdim[1]
#   mu=log(N/sum(exp(en.filter$v+0.5*sigma2true)*xcell*ycell))
#   xcol=seq(xcell/2,plotdim[1],xcell)
#   yrow=seq(ycell/2,plotdim[2],ycell)
#   retry=TRUE
#   r=seq(0,100,1)
#   while(retry){
#     if(nu!=Inf){
#       Y <- GaussRF(x=xcol, y=yrow, model="matern", grid=TRUE,
#                    param=c(mean=0.0, variance=sigma2true, nugget=0.0, scale=alpha,nu=nu))
#     }else{
#       Y <- GaussRF(x=xcol, y=yrow, model="gauss", grid=TRUE,
#                    param=c(mean=0.0, variance=sigma2true, nugget=0.0, scale=alpha/sqrt(0.5)))
#     }
#     Yim <- as.im(list(x = xcol, y = yrow, z = Y))
#     #mu2=log(N/sum(exp(Yim$v)*xcell*ycell))
#     #error=checksimpp(Yim,r,c(sigma2,alpha),nu,mu2)
#     #if(error<2){
#     retry=FALSE
#     #}
#   }
#   
#   Lambda=en.filter
#   Lambda$v=exp(mu+en.filter$v+Yim$v)
#   #simulate inhomogeneous Poisson process with intensity function given by Lambda.
#   X=rpoispp(Lambda)
#   return(X)
# }
# 
# checksimpp=function(Yim,r,par,nu,mu){
#   lambda=Yim
#   lambda$v=exp(mu+Yim$v)
#   X=rpoispp(lambda)
#   re=pcf(X,r=r,correction="translate")
#   re0=pcfmodel(r,par,nu)
#   error=sum((re$trans[-1]-re0[-1])^2)
#   return(error)
# }
# 
# pcfmodel <- function(r, par,nu) {
#   if(nu==Inf){
#     exp(par[1]*exp(-(r/par[2]/sqrt(2))^2))
#   }else
#     exp(Covariance(r,model="matern",param=c(mean=0.0,variance=par[1],nugget=0.0,
#                                             scale=par[2],nu=nu)))
# }
# 
# 
# Sidak <- function(vecP)
#   #
#   # This function corrects a vector of probabilities for multiple testing
#   # using the Bonferroni (1935) and Sidak (1967) corrections.
#   
# {
#   k = length(vecP)
#   
#   vecPB = 0
#   vecPS = 0
#   
#   for(i in 1:k) {
#     bonf = vecP[i]*k
#     if(bonf > 1) bonf=1
#     vecPB = c(vecPB, bonf)
#     vecPS = c(vecPS, (1-(1-vecP[i])^k))
#   }
#   #
#   return(SidakP=vecPS[-1])
# }
# 
# 
# 
# 
# ######################################################################
# calc.pval <- function(ui)
# {
#   # fetch the ui for the observed pattern
#   my.ui <- ui[1]
#   # calculate the number of values greater than the observed
#   my.pval <- (length(ui[ui>=my.ui]))/(length(ui))
#   
#   return(my.pval)
# }
# 
# #######################################################################
# step.ui <- function(Kit, delta.t, npp, c=0.25)
# {
#   Ki.bar <- (sum(Kit) - Kit)/(npp-1)
#   
#   this.ui <- (Kit^c-Ki.bar^c) *
#     (Kit^c-Ki.bar^c) * delta.t
#   
#   
#   ui.step <- vector("numeric", npp)
#   ui.step <- ui.step + this.ui
#   
#   return(ui.step)
# }
# 
# 
# #######################
# 
# formu=function(name,del=NULL){
#   l=length(name)
#   if(l==0)
#     return(as.formula("~1"))
#   
#   for (i in 1:l){
#     if (i==1) fu=paste("~",name[i])
#     else fu=paste(fu,"+",name[i])
#   }
#   if (!is.null(del))
#     fu=paste(fu,del)
#   
#   fu=as.formula(fu)
#   return(fu)
# }
# 
# 
# 
# rMaternppp=function(covr,plotdim,plot=FALSE,
#                     ctrpars=list("nbeta"=6,"Nrange"=c(150,1000),"alpharange"=c(10,60),
#                                  "sigma2range"=c(0.2,2),
#                                  "nurange"=c(0.5,Inf),"missing"=FALSE,
#                                  "nonlinear"=FALSE,"ntry"=30,"delay"=60000)){
#   cln=length(covr)
#   N=round(runif(1,ctrpars$Nrange[1],ctrpars$Nrange[2]),0)
#   alpha=runif(1,ctrpars$alpharange[1],ctrpars$alpharange[2])
#   nu=sample(ctrpars$nurange,1)
#   numCov=ctrpars$nbeta
#   #control tildeZtrue nearly uniform distributed in the range (0,3)
#   ij=sample(1:cln,numCov)
#   beta=runif(numCov,0,2)
#   if(ctrpars$nonlinear){
#     nonij=sample(1:numCov,1)
#   }else{
#     nonij=NA
#   }
#   en.filter=covr[[1]]
#   en.filter$v=covr[[ij[1]]]$v*beta[1]
#   for(i in 2:numCov){
#     if(ctrpars$nonlinear & i==nonij){
#       en.filter$v=en.filter$v+((covr[[ij[i]]]$v)^2)*beta[i]
#     }else{
#       en.filter$v=en.filter$v+covr[[ij[i]]]$v*beta[i]
#     }
#   }
#   en.filter$v=en.filter$v*(2.3)
#   beta=beta*2.3
#   #plot(en.filter)
#   tildeZtrue=var(as.numeric(en.filter$v))
#   sigma2true=runif(1,ctrpars$sigma2range[1],ctrpars$sigma2range[2])
#   
#   parbeta=rep(0,cln)
#   parbeta[ij]=beta
#   partrue=c(tildeZtrue,sigma2true,alpha,round(nu,3),parbeta,round(N,0))
#   X=rLGCP(en.filter,sigma2true,alpha,N,plotdim,nu)
#   itry=0
#   t1=Sys.time()
#   while( (X$n<(N*0.90) | X$n>(N*1.1) ) & itry < ctrpars$ntry ){
#     if(Sys.time()<t1+ctrpars$delay){
#       X=rLGCP(en.filter,sigma2true,alpha,N,plotdim,nu)
#       itry=itry+1
#     }else{
#       itry=ctrpars$ntry+1
#     }
#     
#   }
#   if(itry>ctrpars$ntry)
#     X=NA
#   attr(X,"pars")=partrue
#   if(plot){
#     plot(en.filter)
#     if(!is.na(X))
#       plot(X,add=T)
#   }
#   
#   return(X)
# }
# 
# 
# rMaternPars=function(covr,plotdim,plot=FALSE,
#                      ctrpars=list("nbeta"=6,"Nrange"=c(150,1000),"alpharange"=c(10,60),
#                                   "sigma2range"=c(0.2,2),"nonlinear"=FALSE,
#                                   "nurange"=c(0.5,Inf),"missing"=FALSE,
#                                   "betarange"=c(-2,2))){
#   cln=length(covr)
#   N=round(runif(1,ctrpars$Nrange[1],ctrpars$Nrange[2]),0)
#   alpha=runif(1,ctrpars$alpharange[1],ctrpars$alpharange[2])
#   nu=sample(ctrpars$nurange,1)
#   numCov=ctrpars$nbeta
#   #control tildeZtrue nearly uniform distributed in the range (0,3)
#   ij=sample(1:cln,numCov)
#   beta=runif(numCov,ctrpars$betarange[1],ctrpars$betarange[2])
#   if(ctrpars$nonlinear){
#     nonij=sample(1:numCov,1)
#   }else{
#     nonij=NA
#   }
#   en.filter=covr[[1]]
#   en.filter$v=covr[[ij[1]]]$v*beta[1]
#   for(i in 2:numCov){
#     if(ctrpars$nonlinear & i==nonij){
#       en.filter$v=en.filter$v+((covr[[ij[i]]]$v)^2)*beta[i]
#     }else{
#       en.filter$v=en.filter$v+covr[[ij[i]]]$v*beta[i]
#     }
#   }
#   en.filter$v=en.filter$v
#   beta=beta
#   #plot(en.filter)
#   tildeZtrue=var(as.numeric(en.filter$v))
#   sigma2true=runif(1,ctrpars$sigma2range[1],ctrpars$sigma2range[2])
#   
#   parbeta=rep(0,cln)
#   parbeta[ij]=beta
#   partrue=c(tildeZtrue,sigma2true,alpha,round(nu,3),parbeta,round(N,0))
#   names(partrue)=c("tildeZ","sigma2Y","alpha","nu",paste("beta",1:length(covr),sep=""),"N")
#   attr(partrue,"en.filter")=en.filter
#   if(ctrpars$missing){
#     del=sample(ij,1)
#     attr(partrue,"missing")=del
#   }
#   return(partrue)
# }
# 
# rPattern=function(en.filter,sigma2true,alpha,N,plotdim,nu,
#                   ntry=30, delay=60000,dedge=30,lowl=0.9,upl=1.1,forcechose=FALSE){
#   
#   newwin=owin(c(dedge,plotdim[1]-dedge),c(dedge,plotdim[2]-dedge))
#   tooedge=TRUE
#   itry=0
#   t1=Sys.time()
#   while((itry==0 || X$n<(N*lowl) || X$n>(N*upl) || tooedge ) & itry < ntry ){
#     
#     if(Sys.time()<t1+delay){
#       X=rLGCP(en.filter,sigma2true,alpha,N,plotdim,nu)
#       tooedge=(N/3*2)>(X[newwin])$n
#       itry=itry+1
#     }else{
#       itry=ntry+1
#     }
#     
#   }
#   
#   if(forcechose)
#     return(X)
#   if(itry>ntry)
#     X=NA
#   return(X)
# }
# 
# conf.inter=function(modelpars,covr,plotdim,nsim=19,N,pars){
#   beta=modelpars[7:length(modelpars)]
#   sigma2=modelpars[4]
#   alpha=modelpars[5]
#   nu=modelpars[6]
#   en.filter=covr[[1]]
#   for(i in 1:length(covr)){
#     if(i ==1){
#       en.filter$v=covr[[i]]$v*beta[i]
#     }else{
#       en.filter$v=en.filter$v+covr[[i]]$v*beta[i]
#     }
#   }
#   
#   PVH=unlist(mclapply(1:nsim,shorttry,en.filter=en.filter,sigma2=sigma2,
#                       alpha=alpha,N=N,plotdim=plotdim,nu=nu,data.ppp=data.ppp,
#                       covr=covr[which(beta!=0)],pars=pars,mc.cores=3))
#   #if(nsim==19){
#   #	conf=c(quantile(PVH))
#   #}else{
#   #	conf=sort(PVH)[c(round(nsim*0.025,0),round(nsim*0.975,0))]
#   #}
#   
#   #attr(conf,"PVHs")=PVH
#   return(PVH)
# }
# 
# shorttry=function(x,en.filter,sigma2,alpha,N,plotdim,nu,data.ppp,covr,pars){
#   print(x)
#   retry=TRUE
#   while(retry){
#     data.ppp=rPattern(en.filter,sigma2,alpha,N,plotdim,nu,
#                       ntry=300, delay=60000000,dedge=30,lowl=0.95,upl=1.05,forcechose=TRUE)
#     fittedModel=try(varDecomp.per(data.ppp,covrs=covr,pars=pars,covrsig=FALSE))
#     if(class(fittedModel)!="try-error"){
#       retry=FALSE
#     }
#   }
#   PVH=fittedModel[1]/sum(fittedModel[1:2])
#   return(PVH)
# }
