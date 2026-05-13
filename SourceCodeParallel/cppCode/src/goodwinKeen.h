#ifndef GKVanilla_h
#define GKVanilla_h

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
#include "ODE.h"

using namespace std;

// DEFAULT GOODWIN-KEEN IMPLEMENTATION
template<typename T>
class GoodwinKeenVanilla : virtual public ODE<T> {
public:
	void makeEventTime(const T t, T* parms, T* y, T* x, T h) override {}

	void makeEventVar(const T t, T* parms, T* y, T* x, T h) override {}

	// GOODWIN-KEEN MODEL DEFINITION
	void Func(const T t, const T* y, const T* parms, T* ydot, T* x) override {
		x[0] = 1.0 - y[0] - parms[4] * y[2];
		x[2] = -parms[8] + parms[9]/pow((1.0 - y[1]), 2.0);
		ydot[0] = (x[2] - parms[0]) * y[0];
		x[1] = parms[5] + exp(parms[6] + parms[7] * x[0]);
		ydot[1] = (x[1]/parms[3] - (parms[0] + parms[1] + parms[2])) * y[1];
		ydot[2] = x[1] - x[0] - y[2] * (x[1]/parms[3] - parms[2]);
		}

//	void Func(const T t, const T* y, const T* parms, T* ydot, T* x) override {
//		x[0] = 1.0 - y[0] - parms[4] * y[2];
//		x[2] = -parms[8] + parms[9]/pow((1.0 - y[1]), 2.0);
//		ydot[0] = y[1];
//		x[1] = parms[5] + exp(parms[6] + parms[7] * x[0]);
//		ydot[1] = (1000/3)*(-y[0] + (1 - pow(y[0], 2))*y[1]);
//		ydot[2] = 0;
//		}
};

// GOODWIN-KEEN BUT WITH TWO EXOGENOUS VARIABLES (THIS IS JUST TO ILLUSTRATE THE USE OF INTERMEDIATE VARIABLES)
template<typename T>
class GoodwinKeenExogVar : virtual public ODE<T> {
public:
	void makeEventTime(const T t, T* parms, T* y, T* x, T h) override {}

	void makeEventVar(const T t, T* parms, T* y, T* x, T h) override {}

	// GOODWIN-KEEN MODEL DEFINITION
	void Func(const T t, const T* y, const T* parms, T* ydot, T* x) override {
		x[0] = this->getExogVar(t, 0);
		x[1] = this->getExogVar(t, 1);
		x[2] = 1.0 - y[0] - parms[4] * y[2];
		x[4] = -parms[8] + parms[9]/pow((1.0 - y[1]), 2.0);
		ydot[0] = (x[4] - parms[0]) * y[0];
		x[3] = parms[5] + exp(parms[6] + parms[7] * x[2]);
		ydot[1] = (x[3]/parms[3] - (parms[0] + parms[1] + parms[2])) * y[1];
		ydot[2] = x[3] - x[2] - y[2] * (x[3]/parms[3] - parms[2]);
		}
};

#endif
