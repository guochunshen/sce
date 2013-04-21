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
}

double meanPhylogeneticDistance(NumericVector row1, NumericVector row2,NumericMatrix dis);
double nearestPhylogeneticDistance(NumericVector row1, NumericVector row2, NumericMatrix dis, bool exclude_conspecifics);
