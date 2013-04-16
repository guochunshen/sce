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
#'        a phylogeney should be given. \code{phy$phylosignal} also should be given to indicate how strong the phylogenetic signal in niche or
#'        abundance should be generated. it is the maximum pvalue of the phylogenetic signal in abundance.
#'        
#' @param intra a list represent spatial pattern of intraspecific individuals. it's \code{type} can be "Poisson" and "cluster". if its
#'        \code{type} equals "cluster", a manten clustering process is specified with parameter \code{alpha}, \code{sigma2} and \code{nu}
#'        should be given too.
#' 
#' @param covr a list or an \code{im} object represent map of covariable. if it is a character, the map of covariable will be generated
#'        accroding to this character. current avaiable choise is list(type="sine",scale=1,strength=5). if it is NULL, it means there is effect of covariables on the distribution
#'        of species.
#' @param niche a character controls how does the niche of species generate. current avaiable choise are NULL, "unif", "physignal". the "physignal" 
#'        means the niche of species will generated accroding to the phylogeny of species. in other words, it means there are strong phylogenetic
#'        singal in the ecological niche of species. if it is NULL, it means that there is no effect of habitat filtering and/or competition on the 
#'        distribution of species
#' 
#' @param competition a list to set competition between individuals. if it is NULL, it means there is no competition will happen. otherwise
#'        this parameter should contain at least \code{beta} and \code{r} to define the competition intensity and competition distance between species.
#'        it is also better to give the number of simulation \code{nrep} (default value 5e5) and number of simulation unit \code{nverb} 
#'        (default value 1e5) to report. addtional parameter like \code{verbose} can be configured to disable verb progress report
#'         only interspecific competition was modeled if \code{intra} equals to FALSE.
#'
#'@details
#'it is possible to missing some species if the expected species abundance is quite low.
#'The "physignal" in niche mean there is a high (R2>80) correlation between species niche differences and their phylogenetic correlation
#'
#'@return
#'a scp object
#' 
#'@examples
#'
#'library(testthat)
#'N=1000
#'S=10
#'win=owin(c(0,100),c(0,100))
#'
#'#here, we will test abundance dsitribution, phylogenetic signal in abundance and 
#'# spatial point pattern and spaital phylogenetic community structure
#'
#'#pure random with uniform abundance distribution
#'com=rCom(N,S,win,ab="unif")
#'#the generate species abundance should have no significant difference with the uniform distribtuion
#'expect_true(sum((com$ab-N/S)^2/(N/S))<qchisq(0.95,com$S))
#'
#'
#'#pure random with lognormal abundance distribtuion
#'com=rCom(N,S,win,ab="lognormal")
#'#log transformed abundance should follow normal distribution
#'expect_true(shapiro.test(log(com$ab))$p.value>0.05)
#'
#'
#'#pure random with strong phylogenetic signal abundance distribution
#'com=rCom(N,S,win,ab="physignal",phy=list(br=rexp,phylosignal=0.05))
#'#should be significant phylogenetic signal in abundance
#'expect_true(phylosig(com$phylo,com$ab,test=TRUE)$P<0.1)
#'#the phylogenetic structure should be random at all scale
#'phypvalue=phyMarkCorr(com,cophenetic(com$phylo),nsim=99,rmax=50)$pvalues
#'expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))>0.05)
#'
#'
#'#pure dispersal limitation with logseries abundance distribution
#'com=rCom(N,S,win,ab="logseries",intra=list(type="cluster",sigma2=10,alpha=10,nu=0.5))
#'sp1=subset(com,com$traits$species==names(which(com$ab==max(com$ab))))
#'#there is significant cluster in the spatial distribution of a species
#'expect_true(dclf.test(sp1$com,pcf,nsim=39,verbose=FALSE)$p.value<0.05)
#'
#'
#'#pure habitat filtering with unfirom ahundance distribution and uniform niche distribution
#'com=rCom(N,S,win,ab="unif",covr=list(type="sin",scale=16,strength=5),niche="unif"
#'         ,phy=list(br=runif,phylosignal=0.05))
#'sp1=subset(com,com$traits$species==names(which(com$ab==max(com$ab))))
#'fitmodel=fitCluster(sp1,~covr,sigTest=TRUE)
#'pvalues=as.numeric(attr(fitmodel,"pvalues"))
#'#and have significant habitat association
#'expect_true(pvalues[2]<0.1)
#'#thus phylogenetic structure should be random at all scale
#'comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=99,rmax=50)
#'phypvalue=comphy$pvalues
#'expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))>0.05)
#'plot(comphy)
#'
#'#pure habitat filtering with strong phylogenetic signal in abundance, 
#'#no compeittion, no phylogenetic signal in niche
#'com=rCom(N,S,win,ab="physignal",phy=list(br=runif,phylosignal=0.05),
#'         covr=list(type="sin",scale=16,strength=5),niche="unif")
#'expect_true(phylosig(com$phylo,com$ab,test=TRUE)$P<0.05)
#'#no signal in species niche
#'expect_true(phylosig(com$phylo,com$niche,test=TRUE)$P>0.05)
#'niched=as.matrix(dist(com$niche[match(com$phylo$tip.label,names(com$niche))]))
#'phyd=cophenetic(com$phylo)
#'cor(niched[lower.tri(niched)],phyd[lower.tri(phyd)])
#'#thus phylogenetic structure should be random at all scale
#'comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=99,rmax=50)
#'phypvalue=comphy$pvalues
#'expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))>0.05)
#'plot(comphy)
#'
#'
#'#pure habitat filtering with strong phylogenetic signal in niche
#'com=rCom(N,S,win,ab="unif",phy=list(br=runif,phylosignal=0.05),covr=list(type="sin",scale=16,strength=5),niche="physignal")
#'#phylogenetic signal in niche
#'expect_true(phylosig(com$phylo,com$niche,test=TRUE,nsim=1E3)$P<0.1)
#'#nonrandom phylogenetic pattern at small scale
#'comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=99,rmax=50)
#'phypvalue=comphy$pvalues
#'expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))<0.05)
#'plot(comphy)
#'
#'#pure interspecific competition with unform abundance distribution, and no phylogenetic signal in niche
#'com=rCom(N,S,win,ab="unif",niche="unif",phy=list(br=runif,phylosignal=0.05),
#'         competition=list(beta=0.9,r=10,nrep=5e5,verbose=FALSE,intra=TRUE))
#'pppdata=com$com
#'marks(pppdata)=com$traits$species
#'sp2name=names(which(com$gamma[1,]==min(com$gamma[1,-1])[1])[1])
#'#significant negative association should be existed between species
#'sp2gof=try(dclf.test(pppdata,pcfcross,i="sp1",j=sp2name,nsim=69,r=seq(0,5,length.out=30),verbose=FALSE))
#'if(!inherits(sp2gof,"try-error"))
#'  expect_true(sp2gof$p.value<0.05)
#'#no significant phylogenetic structure is expected at small scale
#'comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=99,rmax=15)
#'phypvalue=comphy$pvalues
#'expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))>0.05)
#'plot(comphy,log="x")
#'
#'
#'#pure competition with unform abundance distribution, 
#'#and signficant phylogenetic signal in niche
#'com=rCom(N,S,win,ab="unif",niche="physignal",phy=list(br=runif,phylosignal=0.05),
#'         competition=list(beta=0.9,r=10,nrep=5e5,verbose=FALSE,intra=TRUE))
#'#phylogenetic signal in niche
#'expect_true(phylosig(com$phylo,com$niche,test=TRUE,nsim=1E3)$P<0.1)
#'#significant phylogenetic structure is expected at small scale
#'comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=199,rmax=30)
#'phypvalue=comphy$pvalues
#'expect_true((1-pbinom(sum(phypvalue<0.05),length(phypvalue),0.05))<0.05)
#'plot(comphy,log="x")
#'
#'
#'#pure competition with strong phylogenetic signal in abundance and niche
#'com=rCom(N,S,win,ab="physignal",niche="physignal",phy=list(br=runif,phylosignal=0.05),
#'         competition=list(beta=0.9,r=5,nrep=5e5,verbose=FALSE,intra=TRUE))
#'expect_true(phylosig(com$phylo,com$ab,test=TRUE)$P<0.05)
#'expect_true(phylosig(com$phylo,com$niche,test=TRUE)$P<0.05)
#'comphy=phyMarkCorr(com,cophenetic(com$phylo),nsim=199,rmax=30)
#'plot(comphy,log="x")
#'


rCom<-function(N,S,win,ab="unif",intra=list(type="Poisson"),phy=NULL,covr=NULL,niche=NULL,competition=NULL, ...){
  
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
    phytree=rcoal(S,tip.label=spname,br=phy$br,...)
  }
  
  #generate the species abundance distribution
  if(ab=="unif"){
    spab=rep(round(N/S),S)
  }else if(ab=="logseries"){
    spab=rlogseries(N,S)
  }else if(ab=="lognormal"){
    spab=rlognormal(N,S)
  }else if(ab=="physignal"){
    spab=fastBM(phytree,a=N/S,sig2=N/S*10,bounds=c(1,N))
    while(phylosignal(spab,phytree)$PIC.variance.P>phy$phylosignal)
      spab=fastBM(phytree,a=N/S,sig2=N/S*10,bounds=c(1,N))
    spab=round(spab)
    spab=spab[match(spname,names(spab))]
  }
  names(spab)=spname
  
  #fix the NA problem in abundance
  spab[is.na(spab)]=1
  
  #generate a ecological relationship if it is needed
  if(!is.null(niche)){
    #if niche is unif, an implicit assumption here is low correlation between niche and phylogeny
    if(niche=="unif"){
      spniche=runif(S)
      if(!is.null(phy)){
        phyd=cophenetic(phytree)
        phyd=phyd[lower.tri(phyd)]
        ror=match(phytree$tip.label,spname)
        niched=as.matrix(dist(spniche[ror]))
        phyniche_cor=cor(phyd,niched[lower.tri(niched)])
        while(phyniche_cor<0 | phyniche_cor>0.001 ){
          spniche=runif(S)
          niched=as.matrix(dist(spniche[ror]))
          phyniche_cor=cor(phyd,niched[lower.tri(niched)])
        }
          
      }
      
    }else if(niche=="physignal"){
#       spniche=fastBM(phytree,a=0.5,bounds=c(0,1))
#       #a significant phylogenetic signal does not mean the correlation between niche distance and phylogenetic distance is high.
#       phyd=cophenetic(phytree)
#       phyd=phyd[lower.tri(phyd)]
#       #ror=match(phytree$tip.label,spname)
#       niched=as.matrix(dist(spniche))
#       phyniche_cor=cor(phyd,niched[lower.tri(niched)])
#       while(phyniche_cor<0.98 ){
#         spniche=fastBM(phytree,a=0.5,bounds=c(0,1))
#         niched=as.matrix(dist(spniche))
#         phyniche_cor=cor(phyd,niched[lower.tri(niched)])
#       }
#       
      spniche=princomp(cophenetic(phytree))$scores[,1]
      
      spniche=scaleRange(spniche,min=0,max=1)
      spniche=spniche[match(spname,names(spniche))]
      
    }
    names(spniche)=spname
  }
  
  #generate a habitat map and intensity map for each species if it is needed
  if(!is.null(covr)){
    if(covr$type=="sin"){
      cl=seq(0,win$xrange[2])
      ncl=length(cl)
      #scaled into range from 0 to 1
      covarable=(sin(cl/covr$scale)+1)/2
      covarable=matrix(rep(covarable,each=ncl),ncol=ncl,nrow=ncl)
    }
    spIntensity=list()
    for(i in 1:S){
      
      tri=exp(-covr$strength*abs(covarable-spniche[i]))
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
    if(niche=="unif"){
      #gamma define the competition strength difference between species based on their niche value
      gamma=as.matrix(dist(spniche))
    }else if(niche=="physignal"){
      gamma=cophenetic(phytree)
      adji=match(spname,rownames(gamma))
      gamma=gamma[adji,]
      gamma=gamma[,adji]
    }
    gamma=round(scaleRange(gamma,0.7,1),6)
    
    rownames(gamma)=colnames(gamma)=spname
    if(!is.null(competition$intra)){
      if(!competition$intra)
        diag(gamma)=1
    }
    #set intraspecific competition is zero here
    
    #radio of interaction
    r=matrix(competition$r,nrow=S,ncol=S)
    rownames(r)=colnames(r)=spname
    spbeta=rep(competition$beta,S)
    names(spbeta)=spname
    mod08 <- list(cif="straussm",par=list(beta=spbeta,gamma=gamma,radii=r),w=win)
    if(!is.null(competition$nrep))
      nrep=competition$nrep
    else
      nrep=5e5
    
    if(!is.null(competition$nverb))
      nverb=competition$nverb
    else
      nverb=1e5
    if(!is.null(competition$verbose)){
      verbose=competition$verbose
    }else
      verbose=TRUE
    if(!verbose){
      nverb=nrep
    }
    #simulated a community just by moving the points.
    com <- rmh(model=mod08,start=list(n.start=spab),control=list(p=1,nrep=nrep,nverb=nverb,fixall=TRUE),verbose=verbose)
    com$marks=data.frame(sp=spname[as.numeric(com$marks)])
    
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
    com$habitat=list(covr=covarable)
    com$niche=spniche
  }
  if(!is.null(competition)){
    com$gamma=gamma
    com$niche=spniche
  }
    
  
  return(com)
  
}