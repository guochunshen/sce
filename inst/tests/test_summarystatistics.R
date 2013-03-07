context("Summary statistic for population distribution")

test_that("the second order statistic",{
  data(testData)
  data.ppp=testData$com
  data.ppm=ppm(data.ppp,~elev+grad,covariates=testData$habitat)
  lambda=predict(data.ppm,locations=data.ppp, type="trend")
  pcf_org=pcf(Kinhom(data.ppp,lambda,r=0:50))
  pcf_adp=pcf_adaptive(data.ppp,lambda,maxt=50,bw=5,adaptive=0.5,kerneltype=1)
})