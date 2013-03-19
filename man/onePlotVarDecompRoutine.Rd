\name{onePlotVarDecompRoutine}
\alias{onePlotVarDecompRoutine}
\title{define a routine function for runing all nessary calcution of variance decompstion for each species in one plot}
\usage{
  onePlotVarDecompRoutine(com, fit_cores = 10,
    gof_cores = 5,
    ctlpars = list(rmax = 25, rmin = 2, bw = 2, sigma2 = 3, alpha = 10, nurange = c(Inf, 0.5), q = 0.5, edgecor = "translate", nsim = 19, r = seq(0, 25, 2), siglevel = 0.05),
    result_rep_id = NULL)
}
\arguments{
  \item{com}{a scp object represent a community or a
  population}

  \item{fit_cores,}{gof_cores number of cores used in model
  fitting and goodness-of-fit test}

  \item{ctlpars}{a list of parameters controled the
  behavior of the model fitting and model diagnoise. Use it
  carefully.}
}
\value{
  a list contains:

  params: parameters of the best fitted cluster model for
  each species. unfitsp: species names that can't be fitted
  under current ctlpars settings. model_performs: pvalues
  of goodness-of_fit test sig_intClu: pvalues of internal
  clustering effect on each species distribution
  sig_habitat: pvalues of habtiat effect on each species
  distribution propVariances: proportions of variance
  explained by each part of the cluster model PVHs:
  proportions of variance explained by habitat PVRs:
  proportions of variance explained by poisson niose(random
  process)
}
\description{
  define a routine function for runing all nessary
  calcution of variance decompstion for each species in one
  plot
}
\details{
  This routine function does four jobs:

  step1: get a list of best fitted cluster models for each
  species step2: evaluate the best fitted model by a
  goodness-of-fit test step3: count how much percentage of
  species significantly affected by habitat and internal
  clustering. step4: calculate variances explained by
  different parts of the model,
}
\examples{
data(testData)
com=testData
com=removeRareSpecies(com,minN=300)

re=onePlotVarDecompRoutine(com)
}

