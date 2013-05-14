#include <Rcpp.h>
#include "myfunction.h"

using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
NumericMatrix commdistInner(const int N, NumericMatrix dis, NumericMatrix x, LogicalMatrix cal_pairs){
  
  NumericMatrix comdist(N,N);
  
  for(int l=0; l<(N-1); l++){
    for(int k=1; k<N; k++){
      if(cal_pairs(k,l) && k>l){
       NumericVector row1=x(k, _);
      NumericVector row2=x(l, _);
      //cout<<N<<" "<<l<<" "<<k<<" ";
      comdist(k,l)=meanPhylogeneticDistance(row1, row2, dis); 
      //cout<<" "<<comdist(k,l)<<endl;
      }
    }
  }
  return comdist;
}


double meanPhylogeneticDistance(NumericVector row1, NumericVector row2, NumericMatrix dis){
          int ncol=row1.size();
          double localsum=0.0;
      
          for(int i=0; i<ncol; i++){
            for(int j=0; j<ncol; j++){
              localsum+=row1[i]*row2[j]*dis(i,j);
            }
          }
         
    return(localsum);
}