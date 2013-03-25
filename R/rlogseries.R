#'
#'  generate a species abundance sample from the logseries distribution given total individual N
#'  and total number of species S
#'  
#'@param N total number of individual
#'@param S total richness   
#'    
#'     
# 

###############################################################################

rlogseries=function(N,S){
	alpha=nlm(lgs,p=1,N=N,S=S,iterlim=1e+4)$estimate
	x=N/(alpha+N)
	ab=numeric()
	for(i in 1:S){
		nmax=(N-sum(ab)-(S-i))
		if(i==S){
			ab[i]=N-sum(ab)
		}else{
			prob=alpha*x^(1:nmax)/(1:nmax)
			prob=prob/sum(prob)
			tp=sample(1:nmax,1,prob=prob)
			while(i<S && tp+sum(ab)>=N)
				tp=sample(1:nmax,1,prob=prob)
			ab[i]=tp
		}
	}
	return(ab)
}

lgs=function(x,N,S){
	if(x<0)
		return(Inf); 
	obs=abs(log(1+N/x)-S/x)
	return(obs)
}


