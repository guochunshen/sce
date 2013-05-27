context("Mark correlation functions")

test_that("the individual based mark correlation function",{
  data(testData)
  #regenerate dbh for each species, every conspecific has the same dbh, thus we can compare our new metric with 
  # the phylogenetic mark correlation function
  spdbh=runif(testData$S)
  testData$traits$dbh=rep(spdbh,testData$ab)
  
  re1=markCorr(com=testData,markName="dbh",r=1:10,testFunName="abdif",nsim=10,h=0.5,exclude_conspecific=TRUE)
  
  spdbh_dist=as.matrix(dist(spdbh))
  rownames(spdbh_dist)=colnames(spdbh_dist)=testData$sp
  
  re2=phyMarkCorr(testData,spdbh_dist,rmax=10,step=1,nsim=10,alpha=0.05,scale=TRUE)
  
  #the observed value should be equal.
  
})