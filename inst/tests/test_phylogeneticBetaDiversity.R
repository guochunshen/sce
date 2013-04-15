context("test on the phylogenetic beta diversity")

test_that("the conventional phylogenetic beta diversity",{
  data(phylocom)
  re1=comdist(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
  #test on the C++ version of the comdist function in the picante package
  re2=comdist_C(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
  #test the result
  expect_equal(re1,re2)
  
  #test the speed advantage
  com=rCom(12000,10,win=owin(c(0,100),c(0,100)),ab="physignal",phy=list(br=rexp,phylosignal=1000))
  com=quadratize(com,20,20)
  sample=table(com$traits$ploti,com$traits$species)
  dist=cophenetic(com$phylo)
  t1=system.time(comdist(sample, dist, abundance.weighted=TRUE))
  t2=system.time(comdist_C(sample, dist, abundance.weighted=TRUE))
  expect_true(as.logical(t1[1]>t2[1]*100))
  
})