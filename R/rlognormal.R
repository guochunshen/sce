#'
#'  generate a lognormal species abundance distribution
#' 
#' @param N total number of individual
#' @param S total number of species
#' @param meanlog the expect mean value of lognormal distribution
#' @param sdlog the expect standard deviation of the lognormal distribution
#' 
#' 
#' 


rlognormal=function(N,S,meanlog=6,sdlog=0.6){
	n=rlnorm(S,meanlog,sdlog)
	n=round(n/sum(n)*N,0)
	if(sum(n)!=N){
		i=which(n==max(n))
		n[i]=N-sum(n)+n[i]
	}
	return(n)
}

