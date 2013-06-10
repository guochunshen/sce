
context("test on the population density estimator")

test_that("test on the sample scheme",{
  #total number of individuals
  N=1000
  #all of the individual distribute in a 1*1 window
  obs=data.frame(x=runif(N),y=runif(N))
  
  #from event to event
  fi=sample(1:N,1)
  focals=obs[fi,]
  for(k in 1:5){
    d1=knq(focals,obs,k=k,q=1,type="etoe")
    d1_obs=sort(sqrt((obs$x[-fi]-focals$x)^2+(obs$y[-fi]-focals$y)^2))[k]
    expect_equal(as.numeric(d1),d1_obs) 
  }
  
  nf=10
  fi=sample(1:N,nf)
  focals=obs[fi,]
  for(k in 1:5){
    d1=knq(focals,obs,k=k,q=1,type="etoe")
    for(nfi in 1:nf){
      d1_obs=sort(sqrt((obs$x[-fi[nfi]]-focals$x[nfi])^2+(obs$y[-fi[nfi]]-focals$y[nfi])^2))[k]
      expect_equal(as.numeric(d1[nfi,]),d1_obs) 
    }
  }
  
  
  
  #from point to event
  focals=data.frame(x=runif(1),y=runif(1))
  for(k in 1:5){
    d1=knq(focals,obs,k=k,q=1,type="ptoe")
    d1_obs=sort(sqrt((obs$x[-fi]-focals$x)^2+(obs$y[-fi]-focals$y)^2))[k]
    expect_equal(as.numeric(d1),d1_obs) 
  }
  
  nf=10
  focals=data.frame(x=runif(nf),y=runif(nf))
  for(k in 1:5){
    d1=knq(focals,obs,k=k,q=1,type="ptoe")
    for(nfi in 1:nf){
      d1_obs=sort(sqrt((obs$x[-fi[nfi]]-focals$x[nfi])^2+(obs$y[-fi[nfi]]-focals$y[nfi])^2))[k]
      expect_equal(as.numeric(d1[nfi,]),d1_obs)
    }
  }
  
  #different sectors
  obs=data.frame(x=c(1,2.1,1.1,2.3,1.5),y=c(1,1,2,2,1.5))
  #from event to event
  focals=obs[5,]
  d2=knq(focals,obs,k=1,q=4,type="etoe")
  d2_obs=sqrt((obs$x[-5]-focals$x)^2+(obs$y[-5]-focals$y)^2)
  expect_equal(sort(d2),sort(d2_obs))
  
  d2=knq(focals,obs,k=1,q=2,type="etoe")
  d2_obs=sqrt((obs$x[-5]-focals$x)^2+(obs$y[-5]-focals$y)^2)
  expect_equal(sort(d2),sort(d2_obs[c(1,3)]))
  
  
  #from point to event
  focals=data.frame(x=1.6,y=1.4)
  d2=knq(focals,obs,k=1,q=4,type="ptoe")
  d2_obs=sqrt((obs$x-focals$x)^2+(obs$y-focals$y)^2)
  expect_equal(d2[1,2],d2_obs[5])
  #test the missing value
  d2=knq(focals,obs,k=2,q=4,type="ptoe")
  
  expect_true(all(is.na(d2[1,-2])))
  expect_equal(d2[1,2],d2_obs[3])
  
  
})


test_that("test on the population density estimator",{
  #total number of individuals
  N=5000
  #all of the individual distribute in a 1*1 window
  obs=data.frame(x=runif(N),y=runif(N))
  
  nf=500
  fi=sample(1:N,nf)
  focals=obs[fi,]
  k=1
  q=1
  type="etoe"
  ds=knq(focals,obs,k=k,q=q,type=type)
  gnonrandomDPDE(ds,dtype=type,k=k,area=1)
  gsimpleDPDE(ds,k=k,area=1)
  
  focals=data.frame(x=runif(nf),y=runif(nf))
  ds=knq(focals,obs,k=k,q=q,type="ptoe")
  gnonrandomDPDE(ds,dtype="ptoe",k=k,area=1)
  gsimpleDPDE(ds,k=k,area=1)
  
  #aggregated distribution
  X <- rThomas(100, 0.1, 30)
  obs=data.frame(x=X$x,y=X$y)
  focals=obs[sample(1:X$n,300),]
  ds=knq(focals,obs,k=k,q=q,type=type)
  gnonrandomDPDE(ds,dtype=type,k=k,area=1)
  
  focals=data.frame(x=runif(nf),y=runif(nf))
  ds=knq(focals,obs,k=k,q=q,type="ptoe")
  gnonrandomDPDE(ds,dtype="ptoe",k=k,area=1)
  gsimpleDPDE(ds,k=k,area=1)
})