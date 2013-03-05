
context("creating and manipulating scp object")

test_that("creating scp object",{
  win=owin(c(0,10),c(0,10))
  sp=sample(c("A","B","C"),50,replace=TRUE)
  x=runif(50,0,50)
  y=runif(50,0,50)
  scpcom=scp(sp,x,y,win)
  
  #check whether it create the same point pattern from input
  expect_equal(scpcom$N,50)
  expect_equal(scpcom$S,3)
  expect_true(all(sort(scpcom$ab)==sort(table(sp))))
  expect_is(scpcom,"scp")
  
  #it should also kept the right relationship between x,y and other point traits. e.g. sp
  expect_is(scpcom$traits,"data.frame")
  expect_equal(scpcom$N,dim(scpcom$traits)[1])
  expect_error(scp(sp,x,y,win,traits=1:10),"traits should be inherits from data.frame")
  expect_error(scp(sp,x,y,win,traits=data.frame(1:10)),"Length of traits not equals to number of individuals")
  
  #the habitat should be a list of im objects
  expect_error(scp(sp,x,y,win,habitat=win),"habitat of the community should be a list of im object")
  vec <- rnorm(1200)
  mat <- matrix(vec, nrow=30, ncol=40)
  whitenoise <- im(mat,xrange=c(0,1),yrange=c(0,1))
  whitenoise2 = im(mat,xrange=c(0,10),yrange=c(0,10))
  expect_error(scp(sp,x,y,win,habitat=list(c(1:3))),"habitat of the community should be a list of im object")
  #it givens warnings if the observed spatial windows is not consistent among different habtiats
  expect_warning(scp(sp,x,y,win,habitat=list(whitenoise,whitenoise2)),"Observed window ranges of habitat are different")
  whitenoise2 = im(mat,xrange=c(0,1),yrange=c(0,1))
  #it gives warnings if the observed spaital windows between individual and habitat are different
  expect_warning(scp(sp,x,y,win,habitat=list(whitenoise,whitenoise2)),"Observed spatial ranges of habitat and individuals are not equal")
  
  #check whether it report equal points warnings
  #maximum resolution of x and y is six
  x[1:2]=c(2,2+1e-7)
  y[1:2]=c(4,4)
  expect_error(scp(sp,x,y,win),"locations of individual are not unique")
  x[1:2]=c(2,2+1e-6)
  y[1:2]=c(4,4)
  expect_equal(scp(sp,x,y,win)$N,50)
  x[1:2]=c(2,2)
  y[1:2]=c(4,4)
  expect_error(scp(sp,x,y,win),"locations of individual are not unique")

  #if we set the forceUnique parameter equals to TRUE, then a uniqued point pattern will be create
  scpcom=scp(sp,x,y,win,forceUnique=TRUE)
  expect_equal(scpcom$N,49)
  
  #check required input parameters
  expect_error(scp(NULL,c(1:4),c(1:4),win),"All of the species, x, y and win parameters are needed")
  expect_error(scp(c(1:4),NULL,c(1:4),win),"All of the species, x, y and win parameters are needed")
  expect_error(scp(c(1:4),c(1:4),NULL,win),"All of the species, x, y and win parameters are needed")
  expect_error(scp(c(1:4),c(1:4),c(1:4),NULL),"All of the species, x, y and win parameters are needed")
  
  #check different length of x,y and traits error
  expect_error(scp(1:3,1:4,1:4,win),"lengthes of species, x and y are not equal")
  expect_error(scp(1:4,1:3,1:4,win),"lengthes of species, x and y are not equal")
  expect_error(scp(1:4,1:4,1:3,win),"lengthes of species, x and y are not equal")
  
  #check type of window object
  expect_error(scp(1:4,1:4,1:4,1:6),"class of win is not owin")
  
  #check infinite or NA values in the x and y
  expect_error(scp(1:4,c(1:3,Inf),1:4,win),"infinite/NA values are not alowed in the coordination of individuals")
  expect_error(scp(1:4,1:4,c(1:3,-Inf),win),"infinite/NA values are not alowed in the coordination of individuals")
  expect_error(scp(1:4,c(NA,1:3),1:4,win),"infinite/NA values are not alowed in the coordination of individuals")
  expect_error(scp(1:4,1:4,c(1:3,NA),win),"infinite/NA values are not alowed in the coordination of individuals")
  
  #check missing NA values in the species name
  expect_warning(scp(c(1:3,NA),1:4,1:4,win),"NA values are contained in species")
  
})

test_that("manipulating scp object",{
  data(testData)
  
  #subset a community
  com1=subset(testData,1:10)
  expect_equal(com1$N,10)
  expect_is(com1$traits$species,"factor")
  expect_equal(nlevels(com1$traits$species),com1$S)
  expect_equal(length(com1$ab),com1$S)

  #select community with specific species abundances
  com2=selSpeciesByAbund(testData,minN=10,maxN=60)
  spab=com2$ab
  expect_true(all(spab>10) & all(spab<60))
  spab2=table(com2$traits$species)
  expect_true(all(sort(spab)==sort(spab2)))

  
})

