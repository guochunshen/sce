#include <Rcpp.h>
#include "myfunction.h"

using namespace Rcpp;


//[[Rcpp::export]]
List phyloBeta_C(const int N, NumericMatrix dis, NumericMatrix x,NumericMatrix qxy, double rmax,
bool exclude_conspecifics=true,int indextype=0){
  
  std::vector<double> spacedists;
  std::vector<double> phydists;
  
  for(int l=0; l<(N-1); l++){
    for(int k=1; k<N; k++){
      if(k>l){
        double pairdist = pow(pow(qxy(k,0)-qxy(l,0),2)+pow(qxy(k,1)-qxy(l,1),2),0.5);
        if(pairdist<=rmax){
          NumericVector row1=x(k, _);
          NumericVector row2=x(l, _);
          double localsum=0.0;
          if(indextype==0){
            localsum=meanPhylogeneticDistance(row1,row2,dis);
          }else if(indextype==1){
            localsum=nearestPhylogeneticDistance(row1,row2,dis,exclude_conspecifics); 
          }
            
          phydists.push_back(localsum);
          spacedists.push_back(pairdist);
          
        }
       
      }
    }
  }
  List results;
  results["spaced"]=spacedists;
  results["phyd"]=phydists;
  return(results);
}
