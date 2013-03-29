/*
 * currently assume conspecific individuals are not interested in the whole analysis
 *
 *  Created on: 2011-12-29
 *      Author: wuping
 */
#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <ctime>

using namespace std;

#include "myfunction.h"

vector<int> spabund2( int* sp, int* ntotal, int* nsp);

double normolizec2(vector<int> spabfocal, vector<int> spabtotal, vector<int> tips, double* phyd);

vector<list<double> > phylocorrsingle2(double* fx, double* fy, int* fsp, double* x, double* y, int* sp,
		double* rmax, double* phyd, int* nfocal, int* ntotal, double* step,
		int* nsp, int* bincout, double* binsum, vector<int> totaltips,
		vector<list<double> > totalsimuvalue, int* nsim, int nsteps, double* normc, double* pk,double* pvalues);

void phylocorr(int* focal, int* nfocal, double* x, double* y, int* ntotal,
		int* sp,int* nsp, double* rmax, double* phyd, double* step, double* pk, int* nsim,
		double* mean, double* upper, double* lower, double* pvalues, double* alpha, int* scale){

	//extract focal individuals information
	double fx[*nfocal];
	double fy[*nfocal];
	int fsp[*nfocal];
	for(int i=0; i< *nfocal; i++){
		fx[i]=x[focal[i]];
		fy[i]=y[focal[i]];
		fsp[i]=sp[focal[i]];
	}

	//abundance of each species
	vector<int> spabfocal = spabund2(fsp,nfocal,nsp);
	vector<int> spabtotal = spabund2(sp,ntotal,nsp);


	//assume tips are start from 0 to nsp, and get nsim number of random permuted tips
	vector<int> tips;
	int totaltipsn=(1+*nsim)*(*nsp);
	vector<int> totaltips;
	for (int i=0; i<*nsp; i++){
		tips.push_back(i);
		totaltips.push_back(i);
	}

	 //calculate normalized c if scale is true
	 double normc[(*nsim+1)];
	 if(*scale==1){
		 normc[0]=normolizec2(spabfocal,spabtotal,tips,phyd);
	 }

	for (int i=1; i< (*nsim+1);i++){
		 //random permute tips
		 random_shuffle(tips.begin(),tips.end());

		 normc[i]=normolizec2(spabfocal,spabtotal,tips,phyd);
		 //cout << normc[i] << " ";
		 for(int j=0; j<*nsp; j++){

			 totaltips.push_back(tips[j]);
		 }
	}
	//for (int i=0; i<totaltips.size();i++)
	//	 cout << totaltips.size();

	//total number of count bins in each simulation
    int nsteps =(int) ceil((*rmax)/(*step));
    int bincout[(*nsim+1)*nsteps];
    double binsum[(*nsim+1)*nsteps];
    vector<list<double> > totalsimuvalue(nsteps);
    //initialize
    for (int i=0; i<(*nsim+1)*nsteps; i++){
    	bincout[i]=0;
    	binsum[i]=0.0;
    }

    totalsimuvalue=phylocorrsingle2(fx, fy, fsp, x, y, sp, rmax, phyd, nfocal, ntotal,
        		step, nsp, bincout, binsum, totaltips, totalsimuvalue,nsim,nsteps,normc,
        		pk,pvalues);


    //number of binsum sets should be recorded
     int nbinset = (int) ceil((*nsim+1)*(*alpha)/2);

    //after the simulation, copy upper and lower limit for retern
    for(int i=0; i<nsteps; i++){
    	totalsimuvalue[i].sort();
    	list<double>::iterator first = totalsimuvalue[i].begin();
    	list<double>::iterator end = totalsimuvalue[i].end();
    	for (int j=0; j<nbinset;j++){
    		first++;
    		end--;
    	}
    	upper[i]=*end;
    	lower[i]=*first;
    	//calculate the both side pvaluves
    	pvalues[i]/=(*nsim+1);
    	if(pvalues[i]>0.5){
    		pvalues[i]=1-pvalues[i];
    	}
    	pvalues[i]*=2;
    }


}


double normolizec2(vector<int> spabfocal, vector<int> spabtotal, vector<int> tips, double* phyd){
	int totalsp = tips.size();
	double phydsum=0.0;
	double totalcounts=0.0;
	for (int i=0; i< spabfocal.size(); i++){
		for (int j=0; j< spabtotal.size(); j++){
			if(i!=j){
				int spjj=tips[j];
				int spii=tips[i];
				int iphy = spjj*(totalsp)+spii;
				//using i's abundance data, and spii's phylogenetic data
				phydsum+=(phyd[iphy]*spabfocal[i]*spabtotal[j]);
				totalcounts+=spabfocal[i]*spabtotal[j];
			}

		}
	}

	return(phydsum/totalcounts);
}

vector<int> spabund2( int* sp, int* ntotal, int*nsp){
	vector<int> spabtotal(*nsp,0);
	for (int j=0; j< *ntotal; j++){
		spabtotal[sp[j]]++;
	}
	return(spabtotal);
}

vector<list<double> > phylocorrsingle2(double* fx, double* fy, int* fsp, double* x, double* y, int* sp,
		double* rmax, double* phyd, int* nfocal, int* ntotal, double* step,
		int* nsp, int* bincout, double* binsum, vector<int> totaltips,
		vector<list<double> > totalsimuvalue, int* nsim, int nsteps, double* normc, double* pk, double* pvalues ){

	double dist=0.0;
	double ibin=0.0;
	int upb=0;
	int lowb=0;
	int iphy=0;

	for (int i = 0; i < *nfocal; i++){
	    for (int j= 0; j < *ntotal; j++){
	    	if( (fsp[i]==sp[j]) ||
	    	    (fx[i]>(x[j]+*rmax)) || (fx[i]<(x[j]-*rmax)) ||
	    	    (fy[i]>(y[j]+*rmax)) || (fy[i]<(y[j]-*rmax)) ||
	    	    (x[j]==fx[i] && y[j]==fy[i] )){
	    	    continue;
	    	}
	    	dist=pow(pow(x[j]-fx[i],2)+pow(y[j]-fy[i],2),0.5);
	    	if(dist > *rmax ){
	    	    continue;
	    	}
	    	ibin=dist / *step;
	    	upb= (int) ceil(ibin)-1;
	    	lowb= (int) floor(ibin);
	    	if(lowb!=0)
	    	    lowb--;
	    	for(int isim=0; isim<(*nsim+1); isim++){
	    		int jumplength=isim*nsteps;
	    		int tipjumps=isim * (*nsp);
	    		 bincout[upb+jumplength]++;
	    		 bincout[lowb+jumplength]++;

	    		 iphy = totaltips[sp[j]+tipjumps]*(*nsp)+totaltips[fsp[i]+tipjumps];
	    		 binsum[upb+jumplength]+=phyd[iphy];
	    		 binsum[lowb+jumplength]+=phyd[iphy];
	    	}
	    }
	}

	for(int i=0; i<(*nsim+1);i++){
		for (int j=0; j< nsteps; j++){
			int locali=j+i*nsteps;
			if(bincout[locali]!=0){
				double localpdi=binsum[locali]/bincout[locali]/normc[i];
				totalsimuvalue[j].push_back(localpdi);
				if(i==0){
					pk[j]=localpdi;
				}else{
					if(pk[j]>localpdi)
						pvalues[j]++;
				}
			}
		}
	}

	return(totalsimuvalue);

}
