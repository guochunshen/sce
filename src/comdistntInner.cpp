#include <Rcpp.h>
#include "myfunction.h"

using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix comdistntInner(const int N, NumericMatrix dis, NumericMatrix x, bool exclude_conspecifics, LogicalMatrix cal_pairs){
  
  NumericMatrix comdist(N,N);
  
  for(int l=0; l<(N-1); l++){
    for(int k=1; k<N; k++){
      if(cal_pairs(k,l) && k>l){
       NumericVector row1=x(k, _);
      NumericVector row2=x(l, _);
     
      comdist(k,l)=nearestPhylogeneticDistance(row1,row2,dis,exclude_conspecifics); 
      }
    }
  }
  return(comdist);
}

double nearestPhylogeneticDistance(NumericVector row1, NumericVector row2, NumericMatrix dis, bool exclude_conspecifics){
   int ncol=row1.size();
      double localsum=0.0;
      
      //from row1 to row2
      for(int i=0; i<ncol; i++){
        if(row1[i]!=0){
          double mindis=0.0;
          bool setvalue=false;
          for(int j=0; j<ncol; j++){
            if(row2[j]!=0){
              //not conspecific
              if(!(exclude_conspecifics && i==j)){
               
                //intial a distance
                if(!setvalue){
                  mindis=dis(i,j);
                  setvalue=true;
                //find the minimum distance
                }else if(dis(i,j)<mindis){
                  mindis=dis(i,j);
                }
                  
              }
            }
          }
        
        localsum+=row1[i]*mindis; 
        }
      }
      
      
      //from row2 to row1
      
      for(int i=0; i<ncol; i++){
        if(row2[i]!=0){
          double mindis=0.0;
          bool setvalue=false;
          for(int j=0; j<ncol; j++){
            if(row1[j]!=0){
              //not conspecific
              if(!(exclude_conspecifics && i==j)){
                //intial a distance
                if(!setvalue){
                  mindis=dis(i,j);
                  setvalue=true;
                //find the minimum distance
                }else if(dis(i,j)<mindis){
                  mindis=dis(i,j);
                }
                  
              }
            }
          }
        localsum+=row2[i]*mindis; 
        }
      }
      return(localsum);
}