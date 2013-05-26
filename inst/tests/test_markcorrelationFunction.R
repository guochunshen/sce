context("Mark correlation functions")

test_that("the individual based mark correlation function",{
  data(testData)
  debug(markCorr)
  re1=markCorr(com=testData,markName="dbh",r=1:10,testFunName="abdif",nsim=10,h=1,exclude_conspecific=FALSE)
  
  #define the test function, m1,m2 is the position of individual in the data set
  data.ppp=testData$com
  marks(data.ppp)=testData$traits$dbh
  isna=is.na(testData$traits$dbh)
  data.ppp=data.ppp[!isna]
  re2=markcorr(data.ppp,f=function(m1,m2) {abs(m1-m2)},r=0:10)
  
  #re2$trans[-1]
  
})