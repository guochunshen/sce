/*
 * mysum.h
 *
 *  Created on: 2011-12-29
 *      Author: wuping
 */
#include <vector>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

extern "C"{

void phylocorr(int* focal, int* nfocal, double* x, double* y, int* ntotal,
		int* sp,int* nsp, double* rmax, double* phyd, double* step, double* pk, int* nsim,
		double* mean, double* upper, double* lower, double* pvalues, double* alpha, int* scale);
void comdist(int* N, int* S, double* dis, double* x, int* cal_pairs, double* comdist);
void comdistnt(int* N, int* S, double* dis, double* x, int* cal_pairs, int* exclude_conspecifics, double* comdist);

}

double meanPhylogeneticDistance(NumericVector row1, NumericVector row2,NumericMatrix dis);
double nearestPhylogeneticDistance(NumericVector row1, NumericVector row2, NumericMatrix dis, bool exclude_conspecifics);
double meanPhylogeneticDistance(int* S,  int* N,int k, int l, double* x, double* dis);
double nearestPhylogeneticDistance(int* S, int* N, int k, int l, double* x, int* exclude_conspecifics,double* dis);