#ifndef GKRK4Vanilla_h
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
#include "goodwinKeen.h"

using namespace std;

template<typename T>
class GoodwinKeenRK4Vanilla: public GoodwinKeenVanilla<T>, public RK4Fixed<T>, virtual public ODE<T> {
public:
	// CONSTRUCTOR
	GoodwinKeenRK4Vanilla(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn) {}	
	//ALTERNATE CONSTRUCTOR (NO DEFAULT VARIABLE)
	GoodwinKeenRK4Vanilla(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn) :
			GoodwinKeenRK4Vanilla{nVIn, nIVIn, tInitIn, tEndIn, ntIn, nullptr} {};
protected:
	exogVar<T>* myExogVar;
};



template<typename T>
class GoodwinKeenRK4EV: public GoodwinKeenExogVar<T>, public RK4Fixed<T>, virtual public ODE<T> {
public:
	// CONSTRUCTOR
	GoodwinKeenRK4EV(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn) {}	
	// Copy constructor
	GoodwinKeenRK4EV(const GoodwinKeenRK4EV<T>& other) : 
		ODE<T>{other.nV, other.nIV, other.tInit, other.tEnd, other.ntIn}, myExogVar(other.myExogVar) {}
protected:
	exogVar<T>* myExogVar;
	T getExogVar(T t, int i) override {
		return myExogVar->getValue(t, i);
	}

};

#endif
