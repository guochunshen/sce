#'
#' define a routine function for runing all nessary calcution of variance decompstion for each species in one plot 
#'
#'@param com  a scp object represent a community or a population
#'@param fit_cores, gof_cores number of cores used in model fitting and goodness-of-fit test
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
#'params: parameters of the best fitted cluster model for each species.
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
#'



onePlotVarDecompRoutine <- function(com,fit_cores=10,gof_cores=5,ctlpars=list(rmax = 25, rmin = 2, 
                          bw = 2, sigma2 = 3, alpha = 10, nurange = c(Inf, 0.5), q = 0.5, 
                            edgecor = "translate",  nsim= 19, r = seq(0, 25, 2), siglevel = 0.05),
                            result_rep_id=NULL){
  
  #run the routine of spatial variance decomposition on each species in each plot
  fitmodel_rep=read_repository(result_rep_id,"./reps")
  if(!inherits(fitmodel_rep,"repository"))
    fitmodel_rep=create_repository(result_rep_id,com$sp,"./reps")
  
  #step1: get a list of best fitted cluster models for each species 
  results=applyGroup(com,com$traits$species,onespModel,ctlpars=ctlpars,verbpro=TRUE,mc.cores=fit_cores,
                          multicore=TRUE,repository=fitmodel_rep)
 
  return(results)
}

#define a convinent function for searching a best fitted cluster model for each species
onespModel=function(com,ctlpars){
  #step1: find the best fitted model
  hnames=names(com$habitat)
  trend=as.formula(paste("~ ", paste(hnames, collapse = "+")))
  fittedmodel=try(fitCluster(com,trend,sigTest=TRUE,ctlpars=ctlpars))
  fittedmodel=try(backwardStep(fittedmodel))
  
  if(inherits(fittedmodel,"try-error"))
    return(NULL)
  
  pvalues=attr(fittedmodel,"pvalues")
  
  #step2 The performance of the best fitted model can be evaluated by a goodness-of-fit test
  model_perform=gofTest(fittedmodel,SFun=Fest,rRange=c(0,ctlpars$rmax),nsim=ctlpars$nsim,r=seq(0,ctlpars$rmax,1),correction="cs")
  
  #step3 calcualte the proportion of variance explained by each part of the model
  propVariance=varDecomp(fittedmodel,r=c(0:ctlpars$rmax),R=4)
  
  #step4 confidence interval of PVH and PVR
  var_conf=envelopeVar(fittedmodel,nsim=ctlpars$nsim,conf_level=0.95,
  r=c(0:ctlpars$rmax),R=4,delta=1,simple=TRUE)
  
  #compose result
  params=as.vector(fittedmodel)
  names(params)=names(fittedmodel)
  
  result=list(params=params,pvalues=pvalues,model_perform=model_perform,propVariance=propVariance,var_conf=var_conf)
  return(result)
}
