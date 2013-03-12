#'
#' define a routine function for runing all nessary calcution of variance decompstion for each species in one plot 
#'
#'@param com  a scp object represent a community or a population
#'@param ctlpars a list of parameters controled the behavior of the model fitting and model diagnoise. Use it carefully.
#'
#'@details
#'This routine function does four jobs:
#'
#'step1: get a list of best fitted cluster models for each species
#'step2: evaluate the best fitted model by a goodness-of-fit test
#'step3: count how much percentage of species significantly affected by habitat and internal clustering.
#'step4: calculate variances explained by different parts of the model, 
#'
#'@return
#'a list contains:
#'
#'fittedmodels: parameters of the best fitted cluster model for each species.
#'unfitsp: species names that can't be fitted under current ctlpars settings.
#'model_performs: pvalues of goodness-of_fit test
#'sig_intClu: pvalues of internal clustering effect on each species distribution
#'sig_habitat: pvalues of habtiat effect on each species distribution
#'propVariances: proportions of variance explained by each part of the cluster model
#'PVHs: proportions of variance explained by habitat
#'PVRs: proportions of variance explained by poisson niose(random process)
#'
#'@examples
#' data(testData)
#' com=testData
#' com=removeRareSpecies(com,minN=300)
#' 
#' re=onePlotVarDecompRoutine(com)
#'



onePlotVarDecompRoutine <- function(com,ctlpars=list(rmax = 25, rmin = 2, 
                          bw = 2, sigma2 = 3, alpha = 10, nurange = c(Inf, 0.5), q = 0.5, 
                            edgecor = "translate",  nsim= 19, r = seq(0, 25, 2), siglevel = 0.05)){
  
  #step1: get a list of best fitted cluster models for each species by applying the above function on each species:
  fittedmodels=applyGroup(com,com$traits$species,onespModel,ctlpars=ctlpars,verbpro=FALSE,mc.cores=20,multicore=TRUE)
  head(fittedmodels)
  #we can still find some species can not be fitted, 
  unfiti=which(unlist(lapply(fittedmodels,function(x) class(x)[1]))=="try-error")
  unfitsp=com$sp[unfiti]
  print(unfitsp)
  if(length(unfitsp)>0)
    fittedmodels=fittedmodels[-unfiti]
  
  #step2: The performance of the best fitted model can be evaluated by a goodness-of-fit test
  model_performs=unlist(mclapply(fittedmodels,FUN=gofTest,SFun=Fest,rRange=c(0,ctlpars$rmax),nsim=ctlpars$nsim,r=seq(0,ctlpars$rmax,1),correction="cs",mc.cores=5))
  #number of models that still not discribed the data well
  sum(model_performs<0.05)
  
  #step3: Based on the best fitted model, we can count how much percentage of species significantly affected by habitat and internal clustering.
  n_sp_fitted=(length(fittedmodels))
  sig_intClu=sum(unlist(lapply(fittedmodels,function(x) {attr(x,"pvalues")[1]<0.05})))/n_sp_fitted
  #define a function for calculate gloabl pvalue of habitat
  global_pvalue=function(pvalues,n,siglevel=0.05){
    return(1-any(pvalues<(1-(1-siglevel)^(1/n))))
  }
  ncovr=length(com$habitat)
  sig_habitat=sum(unlist(lapply(fittedmodels,function(x) {global_pvalue(attr(x,"pvalues")[-1],ncovr)<0.05})))/n_sp_fitted
  
  
  #step4: Now calculate the proportions of variance explained by habitat and internal clustering based on the best fitted models
  if(!is.null(fittedmodels)){
    propVariances=lapply(fittedmodels, varDecomp, r=c(0:ctlpars$rmax),R=4)
    PVHs=unlist(lapply(propVariances,function(x) x[1]))
    PVRs=unlist(lapply(propVariances,function(x) x[2]))
  }
  
  #step5: Finally, thinning and return result
  fittedmodels=lapply(fittedmodels,function(x) {xname=names(x);x=as.vector(x);names(x)=xname;return(x)})
  
  result=list(fittedmodels=fittedmodels,unfitsp=unfitsp,model_performs=model_performs,
              sig_intClu=sig_intClu,sig_habitat=sig_habitat,propVariances=propVariances,
              PVHs=PVHs,PVRs=PVRs)
  return(result)
}

#define a convinent function for searching a best fitted cluster model for each species
onespModel=function(com,ctlpars){
  hnames=names(com$habitat)
  trend=as.formula(paste("~ ", paste(hnames, collapse = "+")))
  fittedmodel=try(fitCluster(com,trend,sigTest=TRUE,ctlpars=ctlpars))
  fittedmodel=try(backwardStep(fittedmodel))
  return(fittedmodel)
}