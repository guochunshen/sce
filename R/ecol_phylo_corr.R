#'
#'generate phylogenetic relationship accroding to the ecological relationship
#'
#'@param ecological_similarity a matrix represent the ecological distance between species
#'@param R2 the expected adjust R square of the correlation between ecological dissimlarity and ecological distance
#'

ecol_phylo_corr=function(ecological_similarity,R2){
  backup=ecological_similarity
  if(is.matrix(ecological_similarity)){
    dimes=dim(ecological_similarity)
    enames=rownames(ecological_similarity)
    diag(ecological_similarity)= NA
    ecological_similarity=as.numeric(ecological_similarity)
    nan=!is.na(ecological_similarity)
    ecological_similarity=ecological_similarity[nan]
  }
  var_es=var(ecological_similarity)
  n=length(ecological_similarity)
  #approximate estimation of the error
  var_err=var_es/R2-var_es
  phylo_dissimilarity=ecological_similarity+rnorm(n,0,var_err)
  phylo_dissimilarity=phylo_dissimilarity-min(phylo_dissimilarity)
  phylo_dissimilarity=phylo_dissimilarity/max(phylo_dissimilarity)
  
  backup[nan]=phylo_dissimilarity
  for(i in 1:(dim(backup)[1])-1){
    for(j in (i+1):(dim(backup)[2])){
      backup[i,j]=backup[j,i]
    }
  }
  
  #dim(phylo_dissimilarity)=dimes
  #rownames(phylo_dissimilarity)=colnames(phylo_dissimilarity)=enames
  return(backup)
}
