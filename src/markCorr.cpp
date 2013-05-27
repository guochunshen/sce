
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

#include "myfunction.h"
//N total number of individuals
//x, y the position of a individual
//mark the mark value of each individual
//r at while place to quantify the mark correlaiton function
//nr length of r
//h torrent bin
//tftype the type of the test function
//nsim number of shuffling used to calcualte confidence interval, 0 means only calculate the observed value
//marki the observed mark index and already permuated index of individual marks, its length equals to N*nsim
//mkcvalues the observed value of the mark correlation function and the simulated value of the mark correlation function under mark shuffling null model



double testfun(double m1, double m2, int* tftype);

void markcorr(int* N, double* x, double* y, double* mark, double* r, int* nr, double* h, int* tftype,
    int* nsim, int* marki, int* exclude_conspecific,int* sp,double* mkcvalues, int* isnormalize, double* tc){
  //define the array to store data
  double tfsum[*nr * (*nsim+1)];
  int pcount[*nr];
  
  double rmax=0;
  //initialize array
  for(int i=0; i< *nr; i++){
    pcount[i]=0;
    for(int j=0; j< (*nsim+1); j++){
        tfsum[i+j * (*nr)]=0;
    }
    
    if(r[i]>rmax){
      rmax=r[i];
    }
  }
  //maximum interpoint distances
  rmax+= *h;
  
  double tcsum=0;
  int tccount=0;
  
  //run over all pairs of points
  for(int i=0; i< *N; i++ ){
    for(int j=i+1; j< *N; j++){
      
      if((*exclude_conspecific)==0 || ((*exclude_conspecific)==1 && sp[i]!=sp[j])){
        double dij=pow(pow(x[i]-x[j],2)+pow(y[i]-y[j],2),0.5);
        //do the calculation only when the interpoint distance is smaller than the maximum r+h
        if(dij<=rmax){
          //one point pair can be assigned into several distance
          vector< int> pis=pointPosition(dij,r,nr,h);
          for (int pij=0; pij<pis.size();pij++ ){
            int pi=pis[pij];
            //if the point is near any r value, just do the point count and testfuntion sum 
            pcount[pi]++;
      
            //sum the test function for observed marks and simulated marks
            for(int k=0; k< (*nsim+1); k++){
              int simui=k* (*nr);
               tfsum[pi+simui]+=testfun(mark[marki[i+simui]],mark[marki[j+simui]], tftype);
            }      
          }
          
        }
        if((*isnormalize)==1){
          tccount++;
          tcsum+=testfun(mark[i],mark[j], tftype);
        }
      }
    }
  }
  //calculate the normalized constant
  if((*isnormalize)==1){
    *tc = tcsum/tccount;
  }else{
    *tc =1;
  }
  //calculate the ratio
  for(int i=0; i< *nr; i++){
    if(pcount[i]!=0){
      
      
      for(int k=0; k< (*nsim+1); k++){
         int simui=k* (*nr);
         mkcvalues[i+simui]=tfsum[i+simui]/pcount[i]/(*tc);
        }
      
      
    }
  }
}

//the test function
double testfun(double m1, double m2, int* tftype){
  double re=0;
  if(*tftype==0){
    // calculate absolute difference between m1 and m2
   re=fabs(m1-m2);
   //just difference from focal mark to neighbor mark
  }else if(*tftype==1){
    re=m1-m2;
    //the sum of marks
  }else if(*tftype==2){
    re=m1+m2;
  }else if(*tftype==3){
    re=fabs(m1-m2)/(m1+m2);
  }
  return re;
}


//it returns the position of a point in the given r vector, and nr otherwise
vector< int> pointPosition(double d, double* r, int* nr, double* h){
  vector< int> pi;
  for(int i=0; i< *nr; i++){
    if(fabs(r[i]-d)<= *h){
      pi.push_back(i);
    }
  }
  return pi;
}