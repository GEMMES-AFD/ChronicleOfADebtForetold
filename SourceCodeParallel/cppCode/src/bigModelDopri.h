#ifndef bigModelDopri_h
#define bigModelDopri_h

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
#include "RK4Fixed.h"
#include "bigModel.h"

using namespace std;

template<typename T>
class bigModelDopri: public bigModel<T>, public dopri<T>, virtual public ODE<T> {
public:
	// CONSTRUCTOR
	bigModelDopri(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn) {}	
	//ALTERNATE CONSTRUCTOR 
	//No exogenous variable
	bigModelDopri(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn) :
		bigModelDopri{nVIn, nIVIn, tInitIn, tEndIn, ntIn, nullptr} {};
	// Constructor with non-default parms for dopri
	bigModelDopri(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn, 
				T hInitIn, T hMinIn, T hMaxIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn), dopri<T>{hInitIn, hMinIn, hMaxIn} {}	
	bigModelDopri(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn, 
				T atolIn, T rtolIn, T facIn, T facMinIn, T facMaxIn, int nStepMaxIn, T hInitIn, T hMinIn, T hMaxIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn), 
			dopri<T>{atolIn, rtolIn, facIn, facMinIn, facMaxIn, nStepMaxIn, hInitIn, hMinIn, hMaxIn} {}	

	// Copy constructor
	bigModelDopri(const bigModelDopri<T>& other) : 
		ODE<T>{other.nV, other.nIV, other.tInit, other.tEnd, other.ntIn}, myExogVar(other.myExogVar), 
		dopri<T>{other.atol, other.rtol, other.fac, other.facMin, other.facMax, other.nStepMax, other.hInit, other.hMax} {}
protected:
	exogVar<T>* myExogVar;
};

template<typename T>
class bigModelRK4: public bigModel<T>, public RK4Fixed<T>, virtual public ODE<T> {
public:
	// CONSTRUCTOR
	bigModelRK4(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn, exogVar<T>* myExogVarIn) :
			ODE<T>{nVIn, nIVIn, tInitIn, tEndIn, ntIn}, myExogVar(myExogVarIn) {}	
	//ALTERNATE CONSTRUCTOR 
	//No exogenous variable
	bigModelRK4(const int nVIn, const int nIVIn, const T tInitIn, const T tEndIn, const int ntIn) :
		bigModelRK4{nVIn, nIVIn, tInitIn, tEndIn, ntIn, nullptr} {};

	// Copy constructor
	bigModelRK4(const bigModelRK4<T>& other) : 
		ODE<T>{other.nV, other.nIV, other.tInit, other.tEnd, other.ntIn}, myExogVar(other.myExogVar) {}
protected:
	exogVar<T>* myExogVar;
};



#endif
