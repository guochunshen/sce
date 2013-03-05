
context("fitting, diagnoising cluster model to species distribution")

test_that("fit cluster model to a species distribution",{
  data(testData)
  sp1=subset(testData,testData$traits$species=="ACALDI")
  #fiting a cluster model without habitat and significant test
  re=fitCluster(sp1)
  expect_is(re,"fm")
  #only have 4 parameters
  expect_equal(length(re),4)
  expect_equal(names(re),c("nu","sigma2","alpha","log(lambda)"))
  expect_equal(attr(re,"pvalues"),NULL)
  
  #fit a cluster model without habitat and do significant test on the internal clustering
  expect_warning(re<-fitCluster(sp1,sigTest=TRUE),"can't test significant of habitat without included it in the model")
  pvalues=attr(re,"pvalues")
  expect_equal(length(pvalues),2)
  #since there is no habitat included in the model, its pvalue should be 1
  expect_equal(as.numeric(pvalues[2]),1)
  #all pvalues should less than or equal to 1
  expect_true(all(pvalues<=1))
  
  #fit a cluster model with habitat
  re<-fitCluster(sp1,~elev+grad)
  #the number of parameters increated to 6 by adding two additional regression coefficients
  expect_equal(length(re),6)
  expect_equal(attr(re,"pvalues"),NULL)
  
  re<-fitCluster(sp1,~elev+grad,sigTest=TRUE)
  #the expected pvalues of aggregated Residual, elev and grad
  expect_equal(length(attr(re,"pvalues")),3)
  expect_equal(names(attr(re,"pvalues")),c("aggreRes","elev","grad"))
  #all pvalues should be less or equals to 1
  expect_true(all(pvalues<=1))
  
  
  #update the fittedmodel automatically
  re<-fitCluster(sp1,~elev+grad,sigTest=TRUE)
  pvalues=attr(re,"pvalues")[-1]
  canimprove=any(pvalues>0.05)
  re_new=updateCluster(re)
  if(canimprove){
    expect_true(length(re)>length(re_new))
    del=which(pvalues==max(pvalues))[1]
    allhabitats=attr(terms(attr(re,"trend")),"term.labels")
    #the habitat with first lagest pvalue will be delete, so its name will not be founded in the results
    expect_true(any(names(re_new)!=allhabitats[del]))
    
  }else{
    expect_true(length(re)==length(re_new))
  }
  #even we didnot do the significant test in the model fitting, it will do it in the update model process
  re2=fitCluster(sp1,~elev+grad,sigTest=FALSE)
  re_new2=updateCluster(re2)
  expect_equal(re_new,re_new2)
  
  
  #if let the program automatically update the model until reached a reasonable model, all habitat pvalues should be signficant
  re_best=backwardStep(re)
  expect_true(all(attr(re_best,"pvalues")[-1]<0.05))
  
})