#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <stdio.h>
#include <string.h>

#include <armadillo>

#if UseRCPP==1
	#include "preprocRCPP_R.h"
#else
	#include "preproc.h"	
#endif

#include "ODE.h"


// TO DO: MAKE A THREAD-SAFE VERSION OF SOLVE TO ALLOW PARALLELIZATION OF CELL MAPPING ALGORITHMS FOR THE RESEARCH OF ZEROES !!!!!!!!!!! 

using namespace std;

template <typename T>
class NR {
public:
	//CONSTRUCTOR
	NR(ODE<T>* aModelIn,
	   T epsilonIn, T tolIn, int nIterMaxIn) :
			myModel(aModelIn),
       			epsilon(epsilonIn),
		        tol(tolIn),
		        nIterMax(nIterMaxIn) {}

	// Numerical approximation of the Jacobian Matrix
	void Jacobian(T t, T* y, T* parms, arma::Mat<T>& out) {
		T yDot0[myModel->getNV()], yDotEpsi[myModel->getNV()], x[myModel->getNIV()];
		myModel->Func(t, y, parms, yDot0, x);
		for (int j=0; j<myModel->getNV(); j++) {
			y[j]+=epsilon;
			myModel->Func(t, y, parms, yDotEpsi, x);
			for(int i=0; i<myModel->getNV(); i++) {
				out(i,j) = (yDotEpsi[i]- yDot0[i])/epsilon; // possible to get rid of this copy if not using arma...
			}
			y[j]-=epsilon;
		}
	}


	// SAME AS ABOVE BUT TAKE DEFAULT VALUES FOR t
	void Jacobian(T* y, T* parms, arma::Mat<T>& out) {
		T t = myModel->getTInit();
		Jacobian(t, y, parms, out);
	}

	// FIND ONE ZERO OF THE SYSTEM USING NEWTON RAPHSON'S ALGORITHM
	void NRMethod(T t, T* y0, T* parms) {
		
		arma::Mat<T> J0(myModel->getNV(), myModel->getNV());
		T yDot0[myModel->getNV()], x[myModel->getNIV()];
		arma::Col<T> Y0(myModel->getNV()), Y1(myModel->getNV()), YDOT0(myModel->getNV());
		int it=0;
		myModel->Func(t, y0, parms, yDot0, x);

		while(norm(yDot0)>tol && it<nIterMax) {
			Jacobian(y0, parms, J0);
			for (int it=0; it<myModel->getNV(); it++) {
				Y0(it) = y0[it];
				YDOT0(it) = -yDot0[it];
			}
			arma::solve(Y1, J0, YDOT0);
			for (int it=0; it<myModel->getNV(); it++) {
				y0[it] = Y1(it) + Y0(it);
			}
			// Update yDot0 at end of loop for stopping condition
			myModel->Func(t, y0, parms, yDot0, x);
			it++;
		}
	};

	// PERFORMS ONE ITERATION OF THE NEWTON RAPHSON ALGORITHM
	void NR1Iter(T t, T* y0, T* parms) {

		T yDot0[myModel->getNV()], x[myModel->getNIV()];
		myModel->Func(t, y0, parms, yDot0, x);
		arma::Mat<T> J0(myModel->getNV(), myModel->getNV());
		Jacobian(y0, parms, J0);
		arma::Col<T> Y0(myModel->getNV()), Y1(myModel->getNV()), YDOT0(myModel->getNV());
		for (int it=0; it<myModel->getNV(); it++) {
			Y0(it) = y0[it];
			YDOT0(it) = -yDot0[it];
		}
		#pragma omp critical 
			{
			arma::solve(Y1, J0, YDOT0);
		}
		if (Y1.n_elem==0) { // if failed to optimize, return arbitrarily low value (to make sur it is sent to sinkcell later)
			for (int it=0; it<myModel->getNV(); it++) y0[it] = -1e250;
			return;
		}
		for (int it=0; it<myModel->getNV(); it++) {
			y0[it] = Y1(it) + Y0(it);
		}
	}
	
	int getNV() {return myModel->getNV();}
	int getTInit() {return myModel->getTInit();}
	int getNRowOut() {return myModel->getNRowOut();}
protected:
	ODE<T>* myModel;
	T norm(T* yDot) {
		T out = 0;
		for (int it=0; it<myModel->getNV(); it++) {
			out+=pow(yDot[it], 2.0);
		}
		return pow(out, 0.5);
	}

	T epsilon;
	T tol;
	int nIterMax;
};



#endif
