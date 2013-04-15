#include <Rcpp.h>

using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix comdistntInner(const int N, NumericMatrix dis, NumericMatrix x, bool exclude_conspecifics){
  
  NumericMatrix comdist(N,N);
  
  for(int l=0; l<(N-1); l++){
    for(int k=1; k<N; k++){
      NumericVector row1=x(k, _);
      NumericVector row2=x(l, _);
      
      int ncol=row1.size();
      double localsum=0.0;
      
      //from i to j
      for(int i=0; i<ncol; i++){
        if(row1[i]!=0){
          double mindis=dis(i,0);
          for(int j=1; j<ncol; j++){
            if(row2[j]!=0 && dis(i,j)<mindis){
              if(!(exclude_conspecifics && i==j)){
                mindis=dis(i,j);
              }
            }
          }
        localsum+=row1[i]*mindis; 
        }
      }
      
      //from j to i
      for(int j=0; j<ncol; j++){
        if(row2[j]!=0){
          double mindis=dis(0,j);
          for(int i=1; i<ncol; i++){
            if(row1[i]!=0 && dis(i,j)<mindis){
             if(!(exclude_conspecifics && i==j)){
                mindis=dis(i,j);
              }
            }
          }
        localsum+=row2[j]*mindis; 
        }
      }
      
      comdist(k,l)=localsum;
    }
  }
  return(comdist);
}