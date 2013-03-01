
#include <Rcpp.h>
#include "stdafx.h"
#include "alglibmisc.h"

using namespace Rcpp;
using namespace alglib;

// [[Rcpp::export]]
List nearest_neighbors(NumericVector xy, int rmax){
  //step1: translate R data formate into alglib data formate
  real_2d_array a;
  int npoints=xy.size()/2;
  a.setcontent(npoints,2,&xy[0]);
  integer_1d_array xytags;
  xytags.setlength(npoints);
  for (int i=0;i<npoints;i++) xytags[i]=i;
  
  
  //step2: build a tagged kdtree
  //define the demontion is 2
  ae_int_t nx = 2;
  ae_int_t ny = 0;
	ae_int_t normtype = 2;
	kdtree kdt;
  kdtreebuildtagged(a, xytags, nx, ny, normtype, kdt);
  
  real_1d_array x;
  ae_int_t k; 
  
  //find the nearest neighbor for each focal point
  integer_1d_array nj = "[]";
  real_1d_array nd = "[]";
  double *a_row=a[0];
   
  for(int i=0; i<npoints;i++){
    x.setcontent(2,a_row);
    k = kdtreequeryrnn(kdt, x, rmax,false);
    std::cout << k<<std::endl;
    kdtreequeryresultstags(kdt, nj);
    kdtreequeryresultsdistances(kdt,nd);
    a_row=a_row+a.getstride();
  }
  
  //transfer result into R formate
  List z = List::create();
  return(z);
}