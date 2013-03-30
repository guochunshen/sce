#' 
#' generate a community with specific properties
#' 
#' @param N expected total number of individuals in the community
#' @param S expected total number of species in the community
#' @param win spaital range of the community
#' @param phy a list contains informations to generate a phylogenetic tree. if it is NULL, it means no phylogenetic information should be simulated
#'        if it is not NULL, there is commonly one parameter is needed to specified a phylogenetic tree: br, the generate function of 
#'        branch length distribution. it can be runif (default) and rexp or other 
#' @param ab a character represent the distribution of species abundance. current values a "unif", "lognormal", "logseries" and "physignal".
#'        note that "physignal" means the species abundance distribution follows phylogenetic relationship between species. in this case,
#'        a phylogeney should be given. \code{phylosignal} also should be given to indicate how strong the phylogenetic signal in niche or
#'        abundance should be generated. 0 means no significant signal. the larger the phylosignal is, the strong the phylogenetic signal is.
#'        
#' @param intra a list represent spatial pattern of intraspecific individuals. it's \code{type} can be "Poisson" and "cluster". if its
#'        \code{type} equals "cluster", a manten clustering process is specified with parameter \code{alpha}, \code{sigma2} and \code{nu}
#'        should be given too.
#' 
#' @param covr a list or an \code{im} object represent map of covariable. if it is a character, the map of covariable will be generated
#'        accroding to this character. current avaiable choise is list(type="sine",scale=1). if it is NULL, it means there is effect of covariables on the distribution
#'        of species.
#' @param niche a character controls how does the niche of species generate. current avaiable choise are NULL, "unif", "physignal". the "physignal" 
#'        means the niche of species will generated accroding to the phylogeny of species. in other words, it means there are strong phylogenetic
#'        singal in the ecological niche of species. if it is NULL, it means that there is no effect of habitat filtering and/or competition on the 
#'        distribution of species
#' 
#' @param competition a list to set competition between individuals. if it is NULL, it means there is no competition will happen. otherwise
#'        this parameter should contain at least \code{beta} and \code{r} to define the competition intensity and competition distance between species.
#'        
#'
#'@details
#'it is possible to missing some species if the expected species abundance is quite low.
#'
#'@return
#'a scp object
#' 
#'@examples
#'
#'N=1000
#'S=10
#'win=owin(c(0,100),c(0,100))
#'
#'#pure random with uniform abundance distribution
#'com=rCom(N,S,win,ab="unif")        
#'com$ab
#'
#'#pure random with lognormal abundance distribtuion
#'com=rCom(N,S,win,ab="lognormal")
#'com$ab
#'
#'#pure random with strong phylogenetic signal abundance distribution
#'com=rCom(N,S,win,ab="physignal",phy=list(br=runif,phylosignal=100))                                                       
#'phylosig(com$phylo,com$ab,test=TRUE)   #pvalue should be very small and below the significant level 
#'
#'#pure dispersal limitation with logseries abundance distribution
#'com=rCom(N,S,win,ab="logseries",intra=list(type="cluster",sigma2=10,alpha=10,nu=0.5))     
#'com$ab                                           
#'sp1=subset(com,com$traits$species==names(which(com$ab==max(com$ab))))
#'plot(sp1)
#'
#'#pure dispersal limitation with strong phylogenetic signal in abundance
#'com=rCom(N,S,win,ab="physignal",intra=list(type="cluster",sigma2=10,alpha=10,nu=0.5),phy=list(br=runif,phylosignal=100))     
#'phylosig(com$phylo,com$ab,test=TRUE)   #pvalue should be very small and below the significant level 
#'sp1=subset(com,com$traits$species==names(which(com$ab==max(com$ab))))
#'plot(sp1)
#'
#'
#'#pure habitat filtering with unfirom ahundance distribution and uniform niche distribution
#'com=rCom(N,S,win,ab="unif",covr=list(type="sin",scale=16),niche="unif")
#'plot(com$habitat)
#'sp1=subset(com,com$traits$species==names(which(com$ab==max(com$ab))))
#'points(sp1$com)
#'
#'
#'#pure habitat filtering with strong phylogenetic signal in abundance, no compeittion
#'com=rCom(N,S,win,ab="physignal",phy=list(br=runif,phylosignal=100),covr=list(type="sin",scale=16),niche="unif")     
#'
#'#pure competition with unform abundance distribution
#'com=rCom(N,S,win,ab="unif",niche="unif",competition=list(beta=0.9,r=5))
#'
#'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
rCom<-function(N,S,win,ab="unif",intra=list(type="Poisson"),phy=NULL,covr=NULL,niche=NULL,competition=NULL){
  
  if(is.null(phy)){
    if(!is.null(niche)){
      if(niche=="physignal")
        stop("can't generate a community with significant phylogenetic signal in niche, but without phylogeney")
    }
    if(ab=="physignal"){
      stop("can't generate a community with significant phylogenetic signal in abundance or niche, but without phylogeney")
    }
  }
  
  if(is.null(niche) & !is.null(competition)){
    stop("can't simulate a competition assembled community without define the niche")
  }
  
  if(is.null(niche) & !is.null(covr))
    stop("can't generate a habitat filtering assembled community without ecological niches")
  
  if(!is.null(covr) & !is.null(competition)){
    stop("Do not support to generate a community with habitat filtering and competition simutanously")
  }
  
  if(intra$type!="Poisson" & !is.null(competition))
    stop("Do not support to generate a intraspecific clustered community with interspecific competition")
  
  #generate species name
  spname=paste("sp",1:S,sep="")
  
  #generate a phylogenetic tree if it is required
  if(!is.null(phy)){
    phytree=rtree(S,tip.label=spname,br=phy$br)
  }
  
  #generate the species abundance distribution
  if(ab=="unif"){
    spab=rep(round(N/S),S)
  }else if(ab=="logseries"){
    spab=rlogseries(N,S)
  }else if(ab=="lognormal"){
    spab=rlognormal(N,S)
  }else if(ab=="physignal"){
    trait1=fastBM(phytree,mu=phy$phylosignal)
    trait1=(trait1-min(trait1))
    trait1=trait1/sum(trait1)
    spab=round(trait1*N)
    spab[spab==0]=1
    spab=spab[match(spname,names(spab))]
  }
  #fix the NA problem in abundance
  spab[is.na(spab)]=1
  
  #generate a ecological relationship if it is needed
  if(!is.null(niche)){
    if(niche=="unif"){
      spniche=runif(S)
    }else if(niche=="physignal"){
      spniche=fastBM(phytree,mu=phy$phylosignal)
      spniche=(spniche-min(spniche))
      spniche=spniche/max(spniche)
      spniche[spniche==0]=1
      spniche=spniche[match(spname,names(spniche))]
    }
  }
  
  #generate a habitat map and intensity map for each species if it is needed
  if(!is.null(covr)){
    if(covr$type=="sin"){
      cl=seq(0,win$xrange[2])
      ncl=length(cl)
      covarable=sin(cl/covr$scale)
      covarable=matrix(rep(covarable,each=ncl),ncol=ncl,nrow=ncl)
    }
    spIntensity=list()
    for(i in 1:S){
      tri=exp(abs(covarable-spniche[i]))
      spIntensity[[i]]=im(tri,xrange=win$xrange,yrange=win$yrange)
    }
    covarable=im(covarable,xrange=win$xrange,yrange=win$yrange)
  }
  
  #generate the point pattern
  
  #case 1: no habitat and no competition, multitype Poisson point process for each species
  if(is.null(covr) & is.null(competition) & intra$type=="Poisson"){
    com=rmpoint(N,win=win,types=spname,ptypes=spab/sum(spab))
    
  #case 2: habitat filtering and no competition, imhomogeneous Poisson process for each species
  }else if(is.null(competition) & !is.null(covr) & intra$type=="Poisson"){
    com=rmpoint(N,f=spIntensity,win=win,types=spname,ptypes=spab/sum(spab))
    
  #case 3: competition and no habitat, a Gibs point process
  }else if(is.null(covr) & !is.null(competition) & intra$type=="Poisson"){
    #gamma define the competition strength difference between species based on their niche value
    gamma=as.matrix(dist(spniche))
    gamma=round(gamma/max(gamma),2)
    gamma=sqrt(gamma+0.5)
    gamma=gamma/max(gamma)
    rownames(gamma)=colnames(gamma)=spname
    
    #radio of interaction
    r=matrix(competition$r,nrow=S,ncol=S)
    mod08 <- list(cif="straussm",par=list(beta=rep(competition$beta,S),gamma=gamma,radii=r),w=win)
    com <- rmh(model=mod08,start=list(n.start=N),
             control=list(ptypes=spab/N,nrep=5e5,nverb=1e5))
    com$marks=data.frame(sp=marks(com))
    #TODO find a way to specify the species name
    
  #case 4: only dispersal limitation, a homogeneous Cox point process fore each species
    #case 5: habitat filtering and dispersal limitation, a imhomogeneous Cox point process for each species  
  }else if(intra$type=="cluster"){
    comx=comy=comsp=c()
    if(is.null(covr))
      en.filter=NULL
    else
      en.filter=spIntensity
    for(i in 1:S){
      spi=rMatern(sigma2=intra$sigma2,alpha=intra$alpha,en.filter=en.filter[[i]],plotdim=c(win$xrange[2],win$yrange[2]),
                  spab[i],nu=intra$nu,xcell=1,ycell=1)
      comx=c(comx,spi$x)
      comy=c(comy,spi$y)
      comsp=c(comsp,rep(spname[i],spi$n))
    }
    com=ppp(x=comx,y=comy,window=win)
    com$marks=data.frame(sp=comsp)
  }
  
  #generating result
  com=as.scp(com)
  if(!is.null(phy)){
    com$phylo=phytree
  }
  if(!is.null(covr)){
    com$habitat=covarable
  }
    
  
  return(com)
  
}