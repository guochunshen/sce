#include <Rcpp.h>

using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix commdistInner(const int N, NumericMatrix dis, NumericMatrix x, LogicalMatrix cal_pairs){
  
  NumericMatrix comdist(N,N);
  
  for(int l=0; l<(N-1); l++){
    for(int k=1; k<N; k++){
      if(cal_pairs(k,l)){
       NumericVector row1=x(k, _);
      NumericVector row2=x(l, _);
      
      int ncol=row1.size();
      double localsum=0.0;
      
      for(int i=0; i<ncol; i++){
        for(int j=0; j<ncol; j++){
          localsum+=row1[i]*row2[j]*dis(i,j);
        }
      }
      
      comdist(k,l)=localsum; 
      }
    }
  }
  return(comdist);
}