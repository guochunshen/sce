context("Mark correlation functions")

test_that("the individual based mark correlation function",{
  
  #generate tha very simple community with 10 individuals belongs to 3 species within a 1 by 1 area
  testCom=scp(species=sample(1:3,10,replace=TRUE),x=runif(10),y=runif(10),win=owin(),traits=data.frame(dbh=runif(10)))
  
  re1=markCorr(com=testCom,markName="dbh",r=seq(0.1,0.5,0.05),testFunName="abdif",nsim=10,h=0.01,exclude_conspecific=FALSE)
  re2=markCorr(com=testCom,markName="dbh",r=seq(0.1,0.5,0.05),testFunName="abdif",nsim=10,h=0.01,exclude_conspecific=TRUE)
  
  #check the normalized term
  ctsum=0
  ctcount=0
  for(i in 1:testCom$N){
    for(j in 1:testCom$N){
      if(i!=j){
        ctcount=ctcount+1
        ctsum=ctsum+abs(testCom$traits$dbh[i]-testCom$traits$dbh[j])
      }
    }
  }
  ctconst=ctsum/ctcount
  expect_equal(ctconst,re1$const)
  
  ctsum2=0
  ctcount2=0
  for(i in 1:testCom$N){
    for(j in 1:testCom$N){
      if(i!=j & testCom$traits$species[i]!=testCom$traits$species[j]){
        ctcount2=ctcount2+1
        ctsum2=ctsum2+abs(testCom$traits$dbh[i]-testCom$traits$dbh[j])
      }
    }
  }
  ctconst2=ctsum2/ctcount2
  expect_equal(ctconst2,re2$const)
  
  #test the observed value between the first value without exclude conspecific
  r1count=0
  r1sum=0
  for(i in 1:testCom$N){
    for(j in 1:i){
      spd=sqrt((testCom$com$x[i]-testCom$com$x[j])^2+(testCom$com$y[i]-testCom$com$y[j])^2)
      if(spd>=(re1$r[1]-0.01) && spd<=(re1$r[1]+0.01) && i!=j){
        r1count=r1count+1
        r1sum=r1sum+abs(testCom$traits$dbh[i]-testCom$traits$dbh[j])
 
      }
    }
  }
  if(r1count==0){
    expect_equal(0,re1$obs[1])
  }else{
    r1value=r1sum/r1count/ctconst
    expect_equal(r1value,re1$obs[1])
    
  }
 
  #test the observed value between the first value
  r2count=0
  r2sum=0
  for(i in 1:testCom$N){
    for(j in 1:i){
      spd=sqrt((testCom$com$x[i]-testCom$com$x[j])^2+(testCom$com$y[i]-testCom$com$y[j])^2)
      if(spd>=(re1$r[4]-0.01) && spd<=(re1$r[4]+0.01) && i!=j){
        r2count=r2count+1
        r2sum=r2sum+abs(testCom$traits$dbh[i]-testCom$traits$dbh[j])
      }
    }
  }
  if(r2count==0){
    expect_equal(0,re1$obs[4])
  }else{
    r2value=r2sum/r2count/ctconst
    expect_equal(r2value,re1$obs[4])
  }
 
  
  #test the observed value between the first value exclude conspecific
  r1count=0
  r1sum=0
  for(i in 1:testCom$N){
    for(j in 1:i){
      spd=sqrt((testCom$com$x[i]-testCom$com$x[j])^2+(testCom$com$y[i]-testCom$com$y[j])^2)
      if(spd>=(re1$r[1]-0.01) && spd<=(re1$r[1]+0.01) && i!=j && testCom$traits$species[i]!= testCom$traits$species[j]){
        r1count=r1count+1
        r1sum=r1sum+abs(testCom$traits$dbh[i]-testCom$traits$dbh[j])
      }
    }
  }
  if(r1count==0){
    expect_equal(0,re2$obs[1])
  }else{
    r1value=r1sum/r1count/ctconst2
    expect_equal(r1value,re2$obs[1])
    
  }
  
  #test the observed value between the first value
  r2count=0
  r2sum=0
  for(i in 1:testCom$N){
    for(j in 1:i){
      spd=sqrt((testCom$com$x[i]-testCom$com$x[j])^2+(testCom$com$y[i]-testCom$com$y[j])^2)
      if(spd>=(re1$r[4]-0.01) && spd<=(re1$r[4]+0.01) && i!=j && testCom$traits$species[i]!= testCom$traits$species[j]){
        r2count=r2count+1
        r2sum=r2sum+abs(testCom$traits$dbh[i]-testCom$traits$dbh[j])
      }
    }
  }
  if(r2count==0){
    expect_equal(0,re2$obs[4])
  }else{
    r2value=r2sum/r2count/ctconst2
    expect_equal(r2value,re2$obs[4])
  }
  
  
})