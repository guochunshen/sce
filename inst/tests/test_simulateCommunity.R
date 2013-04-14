context("simulate communities with different ecological properties and phylogenetic relationship")

test_that("test the previous rCommunity",{
  #let make sure the phylogenetic mark correlation function works correctly here
  com_control=list(S=10,N=1000,ab="lognormal",win=owin(c(0,100),c(0,100)),
                   intra="Poisson",inter="independent",R2=0.9,covr=FALSE,beta=2,r=5,covrscale=16)
  
  #Poisson distribution
  com=rCommunity(com_control)
  #random phylogenetic structure is expected
  phy=phyMarkCorr(as.scp(com$data.ppp),ecol_phylo_corr(com$ecold,0.9),rmax=20,nsim=99)
  
  #dispersal
  com_control$intra="dispersal"
  com=rCommunity(com_control)
  phy=phyMarkCorr(as.scp(com$data.ppp),ecol_phylo_corr(com$ecold,0.9),rmax=20,nsim=99)
  
  #habitat
  com_control$intra="Poisson"
  com_control$covr=TRUE
  com_control$covrscale=5
  com=rCommunity(com_control)
  phy=phyMarkCorr(as.scp(com$data.ppp),ecol_phylo_corr(com$ecold,0.9),rmax=50,nsim=99)
  phy2=phyMarkCorr(as.scp(com$data.ppp),ecol_phylo_corr(com$ecold,0.0001),rmax=50,nsim=99)
  
  
  
  #compeition, relatively time consuming
  com_control$covr=FALSE
  com_control$inter="dependent"
  com=rCommunity(com_control)
  phy=phyMarkCorr(as.scp(com$data.ppp),ecol_phylo_corr(com$ecold,0.9),rmax=20,nsim=99)
  #plot(phy)
  
  #just from the plot(phy), you can see that the phyMarkCorr function works correctly.
})


test_that("test on the rCom",{
  library(phytools)
  N=1000
  S=10
  win=owin(c(0,100),c(0,100))
  
  #here, we will test abundance dsitribution, phylogenetic signal in abundance and 
  # spatial point pattern and spaital phylogenetic community structure
  
  #pure random with uniform abundance distribution
  com=rCom(N,S,win,ab="unif")
  #the generate species abundance should have no significant difference with the uniform distribtuion
  expect_true(sum((com$ab-N/S)^2/(N/S))<qchisq(0.95,com$S))
  
  
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
  com=rCom(N,S,win,ab="unif",covr=list(type="sin",scale=16,strength=5),niche="unif"
           ,phy=list(br=runif,phylosignal=100))
  sp1=subset(com,com$traits$species==names(which(com$ab==max(com$ab))))
  fitmodel=fitCluster(sp1,~covr,sigTest=TRUE)
  pvalues=as.numeric(attr(fitmodel,"pvalues"))
  #and have significant habitat association
  expect_true(pvalues[2]<0.1)
  #thus phylogenetic structure should be random at all scale
  comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=99,rmax=50)
  phypvalue=comphy$pvalues
  expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))>0.05)
  #plot(comphy)
  
  #pure habitat filtering with strong phylogenetic signal in abundance, 
  #no compeittion, no phylogenetic signal in niche
  com=rCom(N,S,win,ab="physignal",phy=list(br=runif,phylosignal=100),
           covr=list(type="sin",scale=16,strength=5),niche="unif")
  expect_true(phylosig(com$phylo,com$ab,test=TRUE)$P<0.05)
  #no signal in species niche
  expect_true(phylosig(com$phylo,com$niche,test=TRUE)$P>0.05)
  niched=as.matrix(dist(com$niche[match(com$phylo$tip.label,names(com$niche))]))
  phyd=cophenetic(com$phylo)
  cor(niched[lower.tri(niched)],phyd[lower.tri(phyd)])
  #thus phylogenetic structure should be random at all scale
  comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=99,rmax=50)
  phypvalue=comphy$pvalues
  expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))>0.05)
  plot(comphy)
  
  
  #pure habitat filtering with strong phylogenetic signal in niche
  com=rCom(N,S,win,ab="unif",phy=list(br=runif,phylosignal=100),covr=list(type="sin",scale=16,strength=5),niche="physignal")
  #phylogenetic signal in niche
  expect_true(phylosig(com$phylo,com$niche,test=TRUE,nsim=1E3)$P<0.1)
  #nonrandom phylogenetic pattern at small scale
  comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=99,rmax=50)
  phypvalue=comphy$pvalues
  expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))<0.05)
  plot(comphy)

  #pure interspecific competition with unform abundance distribution, and no phylogenetic signal in niche
  com=rCom(N,S,win,ab="unif",niche="unif",phy=list(br=runif,phylosignal=100),
           competition=list(beta=0.9,r=10,nrep=5e5,verbose=FALSE,intra=TRUE))
  pppdata=com$com
  marks(pppdata)=com$traits$species
  sp2name=names(which(com$gamma[1,]==min(com$gamma[1,-1])[1])[1])
  #significant negative association should be existed between species
  sp2gof=try(dclf.test(pppdata,pcfcross,i="sp1",j=sp2name,nsim=69,r=seq(0,5,length.out=30),verbose=FALSE))
  if(!inherits(sp2gof,"try-error"))
    expect_true(sp2gof$p.value<0.05)
  #no significant phylogenetic structure is expected at small scale
  comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=99,rmax=15)
  phypvalue=comphy$pvalues
  expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))>0.05)
  #plot(comphy,log="x")
  
  
  #pure competition with unform abundance distribution, 
  #and signficant phylogenetic signal in niche
  com=rCom(N,S,win,ab="unif",niche="physignal",phy=list(br=rexp,phylosignal=100),
           competition=list(beta=0.9,r=5,nrep=5e5,verbose=FALSE,intra=TRUE))
  #phylogenetic signal in niche
  expect_true(phylosig(com$phylo,com$niche,test=TRUE,nsim=1E3)$P<0.1)
  #significant phylogenetic structure is expected at small scale
  comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=199,rmax=20,step=0.5)
  phypvalue=comphy$pvalues
  expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))<0.05)
  par(mfrow=c(1,2))
  plot(comphy,log="x")
  abline(v=2)
  plot(com,col=com$traits$species)
  
  #pure competition with strong phylogenetic signal in abundance and niche
  com=rCom(N,S,win,ab="physignal",niche="physignal",phy=list(br=runif,phylosignal=100),
           competition=list(beta=0.9,r=5,nrep=5e5,verbose=FALSE,intra=TRUE))
  expect_true(phylosig(com$phylo,com$ab,test=TRUE)$P<0.05)
  expect_true(phylosig(com$phylo,com$niche,test=TRUE)$P<0.05)
  #significant phylogenetic structure at small scale
  comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=199,rmax=30)
  phypvalue=comphy$pvalues
  expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))<0.05)
  #plot(comphy,log="x")
  
  
})