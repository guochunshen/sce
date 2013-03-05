#' calculate the species-area relationship
#' 
#' @param currently assume scp is a presence-absence dataframe which contains long, and lat columns. each of other columns are species. each row is a quadrat.
#' 

SAR<-function(scp,area=1:100,nrep=10){
  SAR_approximate(scp,area,nrep)
}

SAR_approximate<-function(scp,area,nrep){
  locationi=which(names(scp)=="long" | names(scp)=="lat")
  locations=scp[,c(locationi)]
  scp=scp[,-c(locationi)]
  #remove some species without any individuals
  occurSum=apply(scp,2,sum,na.rm=TRUE)
  scp=scp[,occurSum!=0]
  #approximatly calculate the region with individuals
  richSum=apply(scp,1,sum,na.rm=TRUE)
  outi=which(richSum==0)
  scp=scp[-outi,]
  locations=locations[-outi,]
  
  #start to sample
  long_all=sort(unique(locations$long))
  lat_all=sort(unique(locations$lat))
  long_n=length(long_all)
  lat_n=length(lat_all)
  
  re_S=matrix(NA,nrow=length(area),ncol=nrep)
  re_area=matrix(NA,nrow=length(area),ncol=nrep)
  max_area=min(c(max(area),long_n,lat_n))
  for(i in 1:max_area){
    ploti=area[i]
    for(j in 1:nrep){
      long_min=sample(1:(long_n-ploti+1),1)
      lat_min=sample(1:(lat_n-ploti+1),1)
      long_max=long_min+ploti-1
      lat_max=lat_min+ploti-1
      sel_plots=locations$long>=long_all[long_min] & locations$long <long_all[long_max] &
        locations$lat>=lat_all[lat_min] & locations$lat < lat_all[lat_max]
      sample_area=sum(sel_plots)
      
      if(sample_area!=0){
        subscp=scp[sel_plots,]
        sample_S=sum(apply(subscp,2,function(x) any(!is.na(x))))
      }else{
        sample_S=0
      }
      re_area[i,j]=sample_area
      re_S[i,j]=sample_S
    }
  }
  mean_S=numeric()
  sd_S=numeric()
  #calcualte mean species richness
  real_max_area=min(max_area^2,dim(scp)[1])
  for(i in 1:real_max_area){
    s_temp=re_S[re_area==i]
    mean_S[i]=mean(s_temp)
    sd_S[i]=sd(s_temp)
  }
  result=data.frame("meanS"=mean_S,"sdS"=sd_S)
  return(result)
}