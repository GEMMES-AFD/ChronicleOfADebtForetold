#ifndef GKEulerVanilla_h
#define GKEulerVanilla_h

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
#include "euler.h"
#include "goodwinKeen.h"

using namespace std;

template<typename T>
class GoodwinKeenEulerVanilla: public GoodwinKeenVanilla<T>, public Euler<T> {
public:
	// CONSTRUCTOR
	GoodwinKeenEulerVanilla(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn) {}	
	//Alternate constructor: no exogenous variable
	GoodwinKeenEulerVanilla(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn) :
			GoodwinKeenEulerVanilla{nVIn, nIVIn, tInitIn, tEndIn, ntIn, nullptr} {};
	
protected:
	exogVar<T>* myExogVar;
};





template<typename T>
class GoodwinKeenEulerEV: public GoodwinKeenExogVar<T>, public Euler<T> {
public:
	// CONSTRUCTOR
	GoodwinKeenEulerEV(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn) {}	
protected:
	exogVar<T>* myExogVar;
	T getExogVar(T t, int i) override {
		return myExogVar->getValue(t, i);
	}
};

#endif
