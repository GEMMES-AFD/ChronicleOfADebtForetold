#ifndef GKRK4iVanilla_h
#define GKRK4Vanilla_h

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
#include "RK4Fixed.h"
#include "ODE.h"

using namespace std;

template<typename T>
class lm : virtual public ODE<T> {
	public:
	void makeEventTime(const T t, T* parms, T* y, T* x) override {} 

	void makeEventVar(const T t, T* parms, T* y, T* x) override {}

	// GOODWIN-KEEN MODEL DEFINITION
	void Func(const T t, const T* y, const T* parms, T* ydot, T* x) override {
		x[0] = parms[0];
		ydot[0] = x[0];
		ydot[1] = parms[1];
		ydot[2] = parms[2];
	}
};

template<typename T>
class lmRK4 : public RK4Fixed<T>, public lm<T>, virtual public ODE<T> {
public:
	// CONSTRUCTOR
	lmRK4(const int nVIn, const int nIVIn, const int ntIn, const T tInitIn, const T tEndIn) :
	      ODE<T>{nVIn, nIVIn}, RK4Fixed<T>{ntIn, tInitIn, tEndIn} {}	
	
protected:
	T getExogVar(T t, int i) override {
		return 0;
	}

};
#endif
