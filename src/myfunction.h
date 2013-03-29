/*
 * mysum.h
 *
 *  Created on: 2011-12-29
 *      Author: wuping
 */
#include <vector>
using namespace std;

extern "C"{

void phylocorr(int* focal, int* nfocal, double* x, double* y, int* ntotal,
		int* sp,int* nsp, double* rmax, double* phyd, double* step, double* pk, int* nsim,
		double* mean, double* upper, double* lower, double* pvalues, double* alpha, int* scale);

}

