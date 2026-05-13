#ifndef BASINBRUTEFORCE_H
#define BASINBRUTEFORCE_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <stdio.h>
#include <string.h>

#if UseRCPP==1
	#include "preprocRCPP_R.h"
#else
	#include "preproc.h"	
#endif

#include "specModel.h"
#include "ODE.h"

using namespace std;

template<typename T>
class basin {
	public:

	// CONSTRUCTOR
	basin(ODE<T>* myModelIn, T* parmsIn, gridModel<T>& myGridIn, convModel<T>& myConvIn) : 
			myModel(myModelIn), parms(parmsIn),
	       		grid(myGridIn.grid), gridSize(myGridIn.gridSize), 
		 	useForConv(myConvIn.useForConv), 
			yeq(myConvIn.yeq), tol(myConvIn.tol), boundsSD(myConvIn.boundsSD) {}

	// COMPUTE DOMAIN OF ATTRACTION FROM GRID grid.
	void computeBasin(int* out, bool useDist=true) {
		if (useDist) {
			computeBasinDist(out);
		} 
		else {
			computeBasinBounds(out);
		}
	}
	void trySolveLastPoint(const T* yInit, const T* parms, T* out) {
		


	}	
	void computeBasinDist(int* out) {
		double yOut[myModel->getNRowOut()];
		#pragma omp parallel for default(shared) private(yOut) if(UseParallel)
		for (int it=0; it<gridSize; it++) {
			try {
				myModel->solveLastPoint(&grid[it*myModel->getNV()], parms, yOut);	
				out[it] = convDist(yOut);
			} catch (...) {
				out[it] = 0; 
			}
		}
	}

	void computeBasinBounds(int* out) {
		double yOut[myModel->getNRowOut()];
		#pragma omp parallel for default(shared) private(yOut) if(UseParallel)
		for (int it=0; it<gridSize; it++) {
			try {
				myModel->solveLastPoint(&grid[it*myModel->getNV()], parms, yOut);	
				out[it] = convBounds(yOut); 
			} catch (...) {
				out[it] = 0; 
			}
		}
	}

	protected:
	ODE<T>* myModel;
	T* parms;

	T* grid;
	int gridSize;	

	const bool* useForConv; 
	const T* yeq;	 
	const T tol;  	
	const T* boundsSD;
	
	// TEST IF TRAJECTORY CONVERGED USING DISTANCE CRITERION
	int convDist(double* y) {
		double distEq=0.0;
		// Test for convergence
		for (int it1 = 0; it1<myModel->getNV(); it1++) {
			if (useForConv[it1]) distEq+= pow(yeq[it1] - y[it1], 2.0);
		}
		distEq = pow(distEq, 0.5);
		if (distEq<tol) {
			return 1;
		}
		return 0;
	}

	// TEST IF TRAJECTORY CONVERGED USING BOUNDARY CRITERION
	int convBounds(double* y) {
		for (int it1 = 0; it1<myModel->getNV(); it1++) {
			if (useForConv[it1] && (isnan(y[it1]) || y[it1]<boundsSD[2*it1] || y[it1]>boundsSD[2*it1+1]))
				return 0;
		}
		return 1;
	}
};
#endif
