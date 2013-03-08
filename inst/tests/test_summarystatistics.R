context("Summary statistic for population distribution")

test_that("the second order statistic",{
  data(testData)

  onesp=subset(testData,testData$traits$species=="MOURMY")
  data.ppp=onesp$com
  data.ppm=ppm(data.ppp,~elev+grad,covariates=testData$habitat)
  lambda=predict(data.ppm,locations=data.ppp, type="trend")
  pcf_org=pcf(Kinhom(data.ppp,lambda,r=0:25))
  pcf_adp=pcf_adaptive(data.ppp,lambda,maxt=25,bw=2,adaptive=0.5,kerneltype=1)
  
})