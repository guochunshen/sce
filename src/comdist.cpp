#include <iostream>

using namespace std;

#include "myfunction.h"
// dis a matrix with S rows and S cols
// x a matrix with N rows and S cols
// cal_pairs a matrix with N rows and N cols
// comdist a matrix with N rows and N cols.

void comdist(int* N, int* S, double* dis, double* x, int* cal_pairs, double* comdist){
  //for each quadrat pair
  for(int l=0; l<(*N-1); l++){
    for(int k=l+1; k< (*N); k++){
      
      int loci=k+l*(*N);
      
      if(cal_pairs[loci] && k>l){
      comdist[loci]=meanPhylogeneticDistance(S, N, k, l, x, dis); 
      }
    }
  }
  
}

double meanPhylogeneticDistance(int* S, int* N, int k, int l, double* x, double* dis){
  double localsum=0.0;
  for(int i=0; i< (*S); i++){
    for(int j=0; j< (*S); j++){
      localsum+=x[k+i*(*N)]*x[l+j*(*N)]*dis[i+j*(*S)];
    }
  }
  return localsum;
}