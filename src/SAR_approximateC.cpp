#include <Rcpp.h> 
#include <stdlib.h>  

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector SAR_approximateC(NumericVector longxy, NumericVector latxy,NumericMatrix scp,
   NumericVector long_all, NumericVector lat_all, IntegerVector area, int nrep, int max_area) {
  int long_n=long_all.size();
  int lat_n=lat_all.size();
  NumericMatrix re_S(area.size(),nrep);
  NumericMatrix re_area(area.size(),nrep);
  
  for(int i=0; i<(max_area-1); i++){
    int ploti = area[i];
    for (int j=0; j<(nrep-1);j++){
      int long_min = rand() % (long_n-ploti+1);
      int lat_min = rand() % (lat_n-ploti+1);
      int long_max=long_min+ploti-1;
      int lat_max=lat_min+ploti-1;
      LogicalVector sel_plots_cond1= longxy>=long_all[long_min] ;
      LogicalVector sel_plots_cond2= longxy <long_all[long_max] ;
      LogicalVector sel_plots_cond3= latxy>=lat_all[lat_min];
      LogicalVector sel_plots_cond4=latxy < lat_all[lat_max];
      for (int k=0; k<sel_plots_cond1.size();k++){
        //TODO implement the unfinished code accroding to SAR_approximate.R
      }
    }
  }
  
  return longxy;
}