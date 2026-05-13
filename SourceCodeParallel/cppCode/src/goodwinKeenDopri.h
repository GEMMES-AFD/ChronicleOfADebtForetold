#ifndef GKDopriVanilla_h
#define GKDopriVanilla_h

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
#include "dopri.h"
#include "goodwinKeen.h"

using namespace std;

template<typename T>
class GoodwinKeenDopriVanilla: public GoodwinKeenVanilla<T>, public dopri<T>, virtual public ODE<T> {
public:
	// CONSTRUCTOR
	GoodwinKeenDopriVanilla(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn) {}	
	//ALTERNATE CONSTRUCTOR 
	//No exogenous variable
	GoodwinKeenDopriVanilla(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn) :
		GoodwinKeenDopriVanilla{nVIn, nIVIn, tInitIn, tEndIn, ntIn, nullptr} {};
	// Constructor with non-default parms for dopri
	GoodwinKeenDopriVanilla(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn, 
				T hInitIn, T hMinIn, T hMaxIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn), dopri<T>{hInitIn, hMinIn, hMaxIn} {}	
	GoodwinKeenDopriVanilla(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn, 
				T atolIn, T rtolIn, T facIn, T facMinIn, T facMaxIn, int nStepMaxIn, T hInitIn, T hMinIn, T hMaxIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn), 
			dopri<T>{atolIn, rtolIn, facIn, facMinIn, facMaxIn, nStepMaxIn, hInitIn, hMinIn, hMaxIn} {}	

	// Copy constructor
	GoodwinKeenDopriVanilla(const GoodwinKeenDopriVanilla<T>& other) : 
		ODE<T>{other.nV, other.nIV, other.tInit, other.tEnd, other.ntIn}, myExogVar(other.myExogVar), 
		dopri<T>{other.atol, other.rtol, other.fac, other.facMin, other.facMax, other.nStepMax, other.hInit, other.hMax} {}
protected:
	exogVar<T>* myExogVar;
};



template<typename T>
class GoodwinKeenDopriEV: public GoodwinKeenExogVar<T>, public dopri<T>, virtual public ODE<T> {
public:
	// CONSTRUCTOR
	GoodwinKeenDopriEV(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn) {}	
	//ALTERNATE CONSTRUCTOR 
	//No exogenous variable
	GoodwinKeenDopriEV(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn) :
		GoodwinKeenDopriEV{nVIn, nIVIn, tInitIn, tEndIn, ntIn, nullptr} {};
	// Constructor with non-default parms for dopri
	GoodwinKeenDopriEV(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn, 
				T hInitIn, T hMinIn, T hMaxIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn), dopri<T>{hInitIn, hMinIn, hMaxIn} {}	
	GoodwinKeenDopriEV(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn, 
				T atolIn, T rtolIn, T facIn, T facMinIn, T facMaxIn, int nStepMaxIn, T hInitIn, T hMinIn, T hMaxIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn), 
			dopri<T>{atolIn, rtolIn, facIn, facMinIn, facMaxIn, nStepMaxIn, hInitIn, hMinIn, hMaxIn} {}	

	// Copy constructor
	GoodwinKeenDopriEV(const GoodwinKeenDopriEV<T>& other) : 
		ODE<T>{other.nV, other.nIV, other.tInit, other.tEnd, other.ntIn}, myExogVar(other.myExogVar), 
		dopri<T>{other.atol, other.rtol, other.fac, other.facMin, other.facMax, other.nStepMax, other.hInit, other.hMax} {}
protected:
	exogVar<T>* myExogVar;
	T getExogVar(T t, int i) override {
		return myExogVar->getValue(t, i);
	}

};

#endif
