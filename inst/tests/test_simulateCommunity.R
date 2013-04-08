context("simulate communities with different ecological properties and phylogenetic relationship")

test_that("test on the rCom",{
  N=1000
  S=10
  win=owin(c(0,100),c(0,100))
  
  #here, we will test abundance dsitribution, phylogenetic signal in abundance and 
  # spatial point pattern and spaital phylogenetic community structure
  
  #pure random with uniform abundance distribution
  com=rCom(N,S,win,ab="unif")
  #the generate species abundance should have no significant difference with the uniform distribtuion
  expect_true(sum((com$ab-100)^2/100)<qchisq(0.95,com$S))
  
  
  #pure random with lognormal abundance distribtuion
  com=rCom(N,S,win,ab="lognormal")
  #log transformed abundance should follow normal distribution
  expect_true(shapiro.test(log(com$ab))$p.value>0.05)
  
  
  #pure random with strong phylogenetic signal abundance distribution
  com=rCom(N,S,win,ab="physignal",phy=list(br=rexp,phylosignal=1000))
  #should be significant phylogenetic signal in abundance
  expect_true(phylosig(com$phylo,com$ab,test=TRUE)$P<0.1)
  #the phylogenetic structure should be random at all scale
  phypvalue=phyMarkCorr(com,cophenetic(com$phylo),nsim=99,rmax=50)$pvalues
  expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))>0.05)
  
  
  #pure dispersal limitation with logseries abundance distribution
  com=rCom(N,S,win,ab="logseries",intra=list(type="cluster",sigma2=10,alpha=10,nu=0.5))
  sp1=subset(com,com$traits$species==names(which(com$ab==max(com$ab))))
  #there is significant cluster in the spatial distribution of a species
  expect_true(dclf.test(sp1$com,pcf,nsim=39,verbose=FALSE)$p.value<0.05)
  
  
  #pure habitat filtering with unfirom ahundance distribution and uniform niche distribution
  com=rCom(N,S,win,ab="unif",covr=list(type="sin",scale=16,strength=5),niche="unif")
  sp1=subset(com,com$traits$species==names(which(com$ab==max(com$ab))))
  fitmodel=fitCluster(sp1,~covr,sigTest=TRUE)
  pvalues=as.numeric(attr(fitmodel,"pvalues"))
  #and have significant habitat association
  expect_true(pvalues[2]<0.1)
  
  
  #pure habitat filtering with strong phylogenetic signal in abundance, no compeittion, no phylogenetic signal in niche
  com=rCom(N,S,win,ab="physignal",phy=list(br=runif,phylosignal=100),covr=list(type="sin",scale=16,strength=5),niche="unif")
  expect_true(phylosig(com$phylo,com$ab,test=TRUE)$P<0.05)
  #no signal in species niche
  expect_true(phylosig(com$phylo,com$niche,test=TRUE)$P>0.05)
  #thus phylogenetic structure should be random at all scale
  comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=99,rmax=50)
  phypvalue=comphy$pvalues
  expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))>0.05)
  
  #here have some error
  #pure habitat filtering with strong phylogenetic signal in niche
  com=rCom(N,S,win,ab="unif",phy=list(br=runif,phylosignal=100),covr=list(type="sin",scale=16,strength=5),niche="physignal")
  #phylogenetic signal in niche
  expect_true(phylosig(com$phylo,com$niche,test=TRUE,nsim=1E4)$P<0.1)
  #nonrandom phylogenetic pattern at small scale
  comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=99,rmax=10)
  phypvalue=comphy$pvalues
  expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))>0.05)
  
  
  #pure interspecific competition with unform abundance distribution
  com=rCom(N,S,win,ab="unif",niche="unif",competition=list(beta=0.9,r=10,nrep=5e5,verbose=FALSE,intra=FALSE))
  pppdata=com$com
  marks(pppdata)=com$traits$species
  sp2name=names(which(com$gamma[1,]==min(com$gamma[1,-1])[1]))
  #significant negative association should be existed between species
  sp2gof=try(dclf.test(pppdata,pcfcross,i="sp1",j=sp2name,nsim=69,r=seq(0,5,length.out=30),verbose=FALSE))
  if(!inherits(sp2gof,"try-error"))
    expect_true(sp2gof$p.value<0.05)
  
  #pure interspecific competition with strong phylogenetic signal in abundance
  com=rCom(N,S,win,ab="physignal",niche="unif",competition=list(beta=0.9,r=10,nrep=5e5,verbose=FALSE,
                              intra=FALSE),phy=list(br=runif,phylosignal=100))
  expect_true(phylosig(com$phylo,com$ab,test=TRUE)$P<0.05)
  
  
  
})