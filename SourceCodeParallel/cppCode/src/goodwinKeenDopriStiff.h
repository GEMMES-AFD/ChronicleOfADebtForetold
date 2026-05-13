#ifndef GKDopriStiff_h
#define GKDopriStiff_h

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
#include "dopriStiffDetect.h"
#include "goodwinKeen.h"

using namespace std;

template<typename T>
class GoodwinKeenDopriStiff: public GoodwinKeenVanilla<T>, public dopriStiff<T>, virtual public ODE<T> {
public:
	// CONSTRUCTOR
	GoodwinKeenDopriStiff(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn) {}
	//ALTERNATE CONSTRUCTOR
	//No exogenous variable
	GoodwinKeenDopriStiff(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn) :
		GoodwinKeenDopriStiff{nVIn, nIVIn, tInitIn, tEndIn, ntIn, nullptr} {};
	// Constructor with non-default parms for dopri
	GoodwinKeenDopriStiff(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn,
				T hInitIn, T hMinIn, T hMaxIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn), dopriStiff<T>{hInitIn, hMinIn, hMaxIn} {}
	GoodwinKeenDopriStiff(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn,
				T atolIn, T rtolIn, T facIn, T facMinIn, T facMaxIn, int nStepMaxIn, T hInitIn, T hMinIn, T hMaxIn, const int nStiffMaxIn, const int nStiffSuccessiveMaxIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn),
			dopriStiff<T>{atolIn, rtolIn, facIn, facMinIn, facMaxIn, nStepMaxIn, hInitIn, hMinIn, hMaxIn, nStiffMaxIn, nStiffSuccessiveMaxIn} {}

	// Copy constructor
	GoodwinKeenDopriStiff(const GoodwinKeenDopriStiff<T>& other) :
		ODE<T>{other.nV, other.nIV, other.tInit, other.tEnd, other.ntIn}, myExogVar(other.myExogVar),
		dopriStiff<T>{other.atol, other.rtol, other.fac, other.facMin, other.facMax, other.nStepMax, other.hInit, other.hMax, other.nStiffMax, other.nStiffSuccessiveMax} {}
protected:
	exogVar<T>* myExogVar;
};



template<typename T>
class GoodwinKeenDopriStiffEV: public GoodwinKeenExogVar<T>, public dopriStiff<T>, virtual public ODE<T> {
public:
	// CONSTRUCTOR
	GoodwinKeenDopriStiffEV(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn) {}
	//ALTERNATE CONSTRUCTOR
	//No exogenous variable
	GoodwinKeenDopriStiffEV(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn) :
		GoodwinKeenDopriStiffEV{nVIn, nIVIn, tInitIn, tEndIn, ntIn, nullptr} {};
	// Constructor with non-default parms for dopri
	GoodwinKeenDopriStiffEV(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn,
				T hInitIn, T hMinIn, T hMaxIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn), dopriStiff<T>{hInitIn, hMinIn, hMaxIn} {}
	GoodwinKeenDopriStiffEV(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn,
				T atolIn, T rtolIn, T facIn, T facMinIn, T facMaxIn, int nStepMaxIn, T hInitIn, T hMinIn, T hMaxIn, const int nStiffMaxIn, const int nStiffSuccessiveMaxIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn),
			dopriStiff<T>{atolIn, rtolIn, facIn, facMinIn, facMaxIn, nStepMaxIn, hInitIn, hMinIn, hMaxIn,nStiffMaxIn, nStiffSuccessiveMaxIn} {}

	// Copy constructor
	GoodwinKeenDopriStiffEV(const GoodwinKeenDopriStiffEV<T>& other) :
		ODE<T>{other.nV, other.nIV, other.tInit, other.tEnd, other.ntIn}, myExogVar(other.myExogVar),
		dopriStiff<T>{other.atol, other.rtol, other.fac, other.facMin, other.facMax, other.nStepMax, other.hInit, other.hMax, other.nStiffMax, other.nStiffSuccessiveMax} {}
protected:
	exogVar<T>* myExogVar;
	T getExogVar(T t, int i) override {
		return myExogVar->getValue(t, i);
	}

};

#endif
