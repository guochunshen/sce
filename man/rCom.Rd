\name{rCom}
\alias{rCom}
\title{generate a community with specific properties}
\usage{
  rCom(N, S, win, ab = "unif",
    intra = list(type = "Poisson"), phy = NULL,
    covr = NULL, niche = NULL, competition = NULL)
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
  should be given. \code{phylosignal} also should be given
  to indicate how strong the phylogenetic signal in niche
  or abundance should be generated. 0 means no significant
  signal. the larger the phylosignal is, the strong the
  phylogenetic signal is.}

  \item{intra}{a list represent spatial pattern of
  intraspecific individuals. it's \code{type} can be
  "Poisson" and "cluster". if its \code{type} equals
  "cluster", a manten clustering process is specified with
  parameter \code{alpha}, \code{sigma2} and \code{nu}
  should be given too.}

  \item{covr}{a list or an \code{im} object represent map
  of covariable. if it is a character, the map of
  covariable will be generated accroding to this character.
  current avaiable choise is list(type="sine",scale=1). if
  it is NULL, it means there is effect of covariables on
  the distribution of species.}

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
  species.}
}
\value{
  a scp object
}
\description{
  generate a community with specific properties
}
\details{
  it is possible to missing some species if the expected
  species abundance is quite low.
}
\examples{
N=1000
S=10
win=owin(c(0,100),c(0,100))

#pure random with uniform abundance distribution
com=rCom(N,S,win,ab="unif")
com$ab

#pure random with lognormal abundance distribtuion
com=rCom(N,S,win,ab="lognormal")
com$ab

#pure random with strong phylogenetic signal abundance distribution
com=rCom(N,S,win,ab="physignal",phy=list(br=runif,phylosignal=100))
phylosig(com$phylo,com$ab,test=TRUE)   #pvalue should be very small and below the significant level

#pure dispersal limitation with logseries abundance distribution
com=rCom(N,S,win,ab="logseries",intra=list(type="cluster",sigma2=10,alpha=10,nu=0.5))
com$ab
sp1=subset(com,com$traits$species==names(which(com$ab==max(com$ab))))
plot(sp1)

#pure dispersal limitation with strong phylogenetic signal in abundance
com=rCom(N,S,win,ab="physignal",intra=list(type="cluster",sigma2=10,alpha=10,nu=0.5),phy=list(br=runif,phylosignal=100))
phylosig(com$phylo,com$ab,test=TRUE)   #pvalue should be very small and below the significant level
sp1=subset(com,com$traits$species==names(which(com$ab==max(com$ab))))
plot(sp1)


#pure habitat filtering with unfirom ahundance distribution and uniform niche distribution
com=rCom(N,S,win,ab="unif",covr=list(type="sin",scale=16),niche="unif")
plot(com$habitat)
sp1=subset(com,com$traits$species==names(which(com$ab==max(com$ab))))
points(sp1$com)


#pure habitat filtering with strong phylogenetic signal in abundance, no compeittion
com=rCom(N,S,win,ab="physignal",phy=list(br=runif,phylosignal=100),covr=list(type="sin",scale=16),niche="unif")
phylosig(com$phylo,com$ab,test=TRUE)


#pure competition with unform abundance distribution
com=rCom(N,S,win,ab="unif",niche="unif",competition=list(beta=0.9,r=5))
com$ab
plot(com$com,col=com$sp)

#pure compeition with strong phylogenetic signal in abundance
com=rCom(N,S,win,ab="physignal",niche="unif",competition=list(beta=0.9,r=5),phy=list(br=runif,phylosignal=100))
phylosig(com$phylo,com$ab,test=TRUE)
}

