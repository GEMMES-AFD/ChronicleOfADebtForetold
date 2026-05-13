#ifndef SPECMODEL_H
#define SPECMODEL_H

#include <iostream>
#include <fstream>
#include <math.h>

#if UseRCPP==1
	#include "preprocRCPP_R.h"
#else 
	#include "preproc.h"
#endif

// STRUCT STORING INPUTS FOR CONVERGENCE DEFINITION
template <typename T>
struct convModel {
	// CONSTRUCTOR
	convModel(const bool useDistIn, 
		  const bool* useForConvIn,
		  const T* yeqIn, 
		  const T tolIn,
		  const T* boundsSDIn) :
			useDist(useDistIn), useForConv(useForConvIn), yeq(yeqIn), tol(tolIn), boundsSD(boundsSDIn) {}
	
	// CONSTRUCTOR FOR DISTANCE CONVERGENCE CRITERIA
	convModel(const bool* useForConvIn, 
		  const T* yeqIn, 
		  const T tolIn) : 
			useDist(true), useForConv(useForConvIn), yeq(yeqIn), tol(tolIn), boundsSD(nullptr) {}
	// CONSTRUCTOR FOR BOUNDS CONVERGENCE CRITERIA
	convModel(const bool* useForConvIn,
		  const T* boundsSDIn) :
			useDist(false), useForConv(useForConvIn), yeq(nullptr), tol(0), boundsSD(boundsSDIn) {}

	
	const bool useDist;	 // Which convergence criteria to use (distance minimization or convergence domain)
	const bool* useForConv;  // Which variables are considered to test convergence criteria (used for both criteria)
	const T* yeq;	 // Equilibrium position of the model
	const T tol;  	 // Tolerance for convergence 
	const T* boundsSD;	 // Lower & upper bounds of convergence domain
};


// STRUCT DEFINING THE GRID TO EXPLORE
template <typename T>
struct gridModel {
	T* grid;
	int gridSize;	
};

// STRUCT STORING INPUTS FOR SENSITIVITY ANALYSIS template <typename T>
template<typename T>
struct parmsSA {
	const int nParmsSetSA;  // Number of sets of parameters to explore
	const int nParmsSA;	// Number of different parameters per set
	const int* parmsPosSA;  // position of parameters on which we perform the sensitivity analysis within the parms vector
	T* allParmsSA;     // Set of parameters for sensitivity analysis
};

#endif
