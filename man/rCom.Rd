\name{rCom}
\alias{rCom}
\title{generate a community with specific properties}
\usage{
  rCom(N, S, win, ab = "unif",
    intra = list(type = "Poisson"), phy = NULL,
    covr = NULL, niche = NULL, competition = NULL, ...)
}
\arguments{
  \item{N}{expected total number of individuals in the
  community}

  \item{S}{expected total number of species in the
  community}

  \item{win}{spaital range of the community}

  \item{phy}{a list contains informations to generate a
  phylogenetic tree. if it is NULL, it means no
  phylogenetic information should be simulated if it is not
  NULL, there is commonly one parameter is needed to
  specified a phylogenetic tree: br, the generate function
  of branch length distribution. it can be runif (default)
  and rexp or other}

  \item{ab}{a character represent the distribution of
  species abundance. current values a "unif", "lognormal",
  "logseries" and "physignal".  note that "physignal" means
  the species abundance distribution follows phylogenetic
  relationship between species. in this case, a phylogeney
  should be given. \code{physigPvalue} also should be given
  to indicate how strong the phylogenetic signal in niche
  or abundance should be generated.}

  \item{intra}{a list represent spatial pattern of
  intraspecific individuals. it's \code{type} can be
  "Poisson" and "cluster". if its \code{type} equals
  "cluster", a manten clustering process is specified with
  parameter \code{alpha}, \code{sigma2} and \code{nu}
  should be given too.}

  \item{covr}{a list or an \code{im} object represent map
  of covariable. if it is a character, the map of
  covariable will be generated accroding to this character.
  current avaiable choise is
  list(type="sine",scale=1,strength=5). if it is NULL, it
  means there is effect of covariables on the distribution
  of species.}

  \item{niche}{a character controls how does the niche of
  species generate. current avaiable choise are NULL,
  "unif", "physignal". the "physignal" means the niche of
  species will generated accroding to the phylogeny of
  species. in other words, it means there are strong
  phylogenetic singal in the ecological niche of species.
  if it is NULL, it means that there is no effect of
  habitat filtering and/or competition on the distribution
  of species}

  \item{competition}{a list to set competition between
  individuals. if it is NULL, it means there is no
  competition will happen. otherwise this parameter should
  contain at least \code{beta} and \code{r} to define the
  competition intensity and competition distance between
  species.  it is also better to give the number of
  simulation \code{nrep} (default value 5e5) and number of
  simulation unit \code{nverb} (default value 1e5) to
  report. addtional parameter like \code{verbose} can be
  configured to disable verb progress report only
  interspecific competition was modeled if \code{intra}
  equals to FALSE.}
}
\value{
  a scp object
}
\description{
  generate a community with specific properties
}
\details{
  it is possible to missing some species if the expected
  species abundance is quite low. The "physignal" in niche
  mean there is a high (R2>80) correlation between species
  niche differences and their phylogenetic correlation
}
\examples{
library(testthat)
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
plot(comphy)

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
plot(comphy,log="x")


#pure competition with unform abundance distribution,
#and signficant phylogenetic signal in niche
com=rCom(N,S,win,ab="unif",niche="physignal",phy=list(br=runif,phylosignal=100),
        competition=list(beta=0.9,r=10,nrep=5e5,verbose=FALSE,intra=TRUE))
#phylogenetic signal in niche
expect_true(phylosig(com$phylo,com$niche,test=TRUE,nsim=1E3)$P<0.1)
#significant phylogenetic structure is expected at small scale
comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=199,rmax=30)
phypvalue=comphy$pvalues
expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))<0.05)
plot(comphy,log="x")


#pure competition with strong phylogenetic signal in abundance and niche
com=rCom(N,S,win,ab="physignal",niche="physignal",phy=list(br=runif,phylosignal=100),
        competition=list(beta=0.9,r=5,nrep=5e5,verbose=FALSE,intra=TRUE))
expect_true(phylosig(com$phylo,com$ab,test=TRUE)$P<0.05)
expect_true(phylosig(com$phylo,com$niche,test=TRUE)$P<0.05)
comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=199,rmax=30)
plot(comphy,log="x")
}

