#include <iostream>

using namespace std;

#include "myfunction.h"
// dis a matrix with S rows and S cols
// x a matrix with N rows and S cols
// cal_pairs a matrix with N rows and N cols
// comdist a matrix with N rows and N cols.

void comdistnt(int* N, int* S, double* dis, double* x, int* cal_pairs, int* exclude_conspecifics, double* comdist){
  //for each quadrat pair
  for(int l=0; l<(*N-1); l++){
    for(int k=l+1; k< (*N); k++){
      
      int loci=k+l*(*N);
      
      if(cal_pairs[loci] && k>l){
        //from row k to row l
      comdist[loci]=nearestPhylogeneticDistance(S, N, k, l, x, exclude_conspecifics, dis); 
      //from row l to row k
      comdist[loci]+=nearestPhylogeneticDistance(S, N, l, k, x, exclude_conspecifics, dis); 
      }
    }
  }
  
}

double nearestPhylogeneticDistance(int* S, int* N, int k, int l, double* x, int* exclude_conspecifics,double* dis){
  double localsum=0.0;
  for(int i=0; i< (*S); i++){
    int r1 = k+i*(*N);
    //if this quadrat contains this species
    if(x[r1]!=0){
      //find the nearest nearest neighborhood distance
      bool setvalue=false;
      double minivalue=0;
      for(int j=0; j< (*S); j++){
        int r2=l+j*(*N);
        //if the second quadrat contains jth spcies
        if(x[r2]!=0){
          if(i!=j || (i==j && (*exclude_conspecifics!=1))){
            if(!setvalue){
              minivalue=dis[i+j*(*S)];
              setvalue=true;
            }else if(minivalue>dis[i+j*(*S)]){
              minivalue=dis[i+j*(*S)];
            }
          }
        }
      
      }
      localsum+=x[r1]*minivalue;
    }
  
  }
  return localsum;
}