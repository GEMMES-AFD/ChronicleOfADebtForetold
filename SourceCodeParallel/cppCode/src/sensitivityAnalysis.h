#ifndef SA_h
#define SA_h

#if UseRCPP==1
	#include "preprocRCPP_R.h"
#else 
	#include "preproc.h"
#endif

#include "specModel.h"
#include "ODE.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <stdio.h>
#include <string.h>

using namespace std;

template<typename T>
class SA {
public:
	// CONSTRUCTOR 
	SA(ODE<T>* aModelIn, T* parmsInitIn, int parmsLengthIn,
	      int nParmsSetSAIn, int nParmsSAIn, int* parmsPosSAIn, T* allParmsSAIn ) : 
	       		myModel(aModelIn), parmsInit(parmsInitIn), parmsLength(parmsLengthIn),
			nParmsSetSA(nParmsSetSAIn), nParmsSA(nParmsSAIn), parmsPosSA(parmsPosSAIn), allParmsSA(allParmsSAIn) {}	
	
	// Runs model for all parms values in allParmsSA and return only last point of trajectory
	void sensitivityAnalysisLastPoint(T* y0, T* out) {
		// init parms vector (hard copy)
		T parms[parmsLength];
		for (int it=0; it<parmsLength; it++) {
			parms[it] = parmsInit[it];
		}
		#pragma omp parallel for firstprivate(parms) if(UseParallel) 
		for (int it=0; it<nParmsSetSA; it++) {
			for (int it2=0; it2<nParmsSA; it2++) {
				parms[parmsPosSA[it2]] = allParmsSA[it*nParmsSA + it2];
			}
			
			myModel->solveLastPoint(y0, parms, &out[it*myModel->getNRowOut()]);
		}
	}

	void sensitivityAnalysisFullTrajectory(T* y0, T* out) {
		
		// init parms vector (hard copy)
		T parms[parmsLength];
		for (int it=0; it<parmsLength; it++) {
			parms[it] = parmsInit[it];
		}

		#pragma omp parallel for firstprivate(parms) if(UseParallel)
		for (int it=0; it<nParmsSetSA; it++) {
			for (int it2=0; it2<nParmsSA; it2++) {
				parms[parmsPosSA[it2]] = allParmsSA[it*nParmsSA + it2];
			}
			
			myModel->solve(y0, parms, &out[it*myModel->getNRowOut()*myModel->getNt()]);
		}
	}

protected:
	ODE<T>* myModel;
	T* parmsInit;
	int parmsLength;
	int nParmsSetSA;
	int nParmsSA;
	int* parmsPosSA;
	T* allParmsSA;
};

#endif
