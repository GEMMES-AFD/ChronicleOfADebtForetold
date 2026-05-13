#ifndef FIXEDOBJFUNC_H
#define FIXEDOBJFUNC_H


#include <armadillo>

#if UseRCPP==1
	#include "preprocRCPP_R.h"
#else
	#include "preproc.h"	
#endif

#include "minDist.h"


// CLASS TO COMPUTE DISTANCE (objective function for calibration) BUT USING ONLY SELECTED PARAMETERS
// DEFAULT VALUES SET IN parmsFullInit ARE USED FOR NON-SELECTED PARAMETERS

template<typename T, class objFunc>
class fixedObjFunc: public objFunc {
public: 
	using objFunc::Evaluate;
	// CONSTRUCTOR
	fixedObjFunc(objFunc& anObjFuncIn, T* parmsFullInitIn, bool* parmsFixedIn) : 
		objFunc{anObjFuncIn}, 
		parmsFullInit(parmsFullInitIn), parmsFixed(parmsFixedIn) {}

	// COPY CONSTRUCTOR (not needed)
	//fixedObjFunc(const fixedObjFunc& anObjFuncIn) :
		//objFunc{anObjFuncIn}, 
	       // parmsFullInit(anObjFuncIn.parmsFullInit), parmsFixed(anObjFuncIn.parmsFixed) {}
	

	// BUILD FULL VECTOR PARAMETERS FROM PARTIAL PARMS VECTOR (DEFINING ONLY FIXED ELEMENTS)	
	void completeParms(T* parmsFull, arma::Col<T> parmsPart) {
		unsigned int compt=0;
		for (unsigned int it=0; it<this->parmsLength; it++) {
			if (parmsFixed[it]==true) {
				parmsFull[it] = parmsPart[compt];
				compt++;
			} else {
				parmsFull[it] = parmsFullInit[it];
			}
		}

	}

	T Evaluate(arma::Col<T> armaParmsPart) {
		//T parmsFull[this->parmsLength];
		T parmsFull[this->parmsLength];
		completeParms(parmsFull, armaParmsPart);
		T out = this->Evaluate(&parmsFull[0]);
		return out;
	}
	
protected:	
	T* parmsFullInit;
	bool* parmsFixed;

};

#endif
