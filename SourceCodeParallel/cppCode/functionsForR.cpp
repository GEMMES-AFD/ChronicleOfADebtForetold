#define UseRCPP 1 		// true if the code is called from R(via RCPP), false if it is directly called from c++


#include "src/preprocRCPP_R.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <Rcpp.h>

#if UseParallel!=0
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
#endif

#define ARMA_WARN_LEVEL 0
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins("cpp14")]]


#include "src/exogenousVariables.h"
#include "src/modelR.h"
#include "src/solverModelR.h"
#include "src/ODE.h"
#include "src/euler.h"
#include "src/sensitivityAnalysis.h"
#include "src/basinBruteForce.h"
#include "src/CMAES.h"
#include "src/minDist.h"
#include "src/fixedObjFunc.h"
#include "src/HCMSubdivision.h"
#include "src/cellMapSubdivisionNR.h"
#include "src/SCM.h"
#include "src/cellMapNR.h"

using namespace std;


//For now, cubic interpolation is imposed by default (in R)
// [[Rcpp::export]]
vector<double> EulerForR(int nt,                                 // number of points to compute for trajectory
		         double tInit,				 // time at beginning of period to simulate
			 double tEnd,				 // time at end of period to simulate
			 int nV,                                 // number of variables of the system
			 int nIV,                                // number of exogenous variables
		         vector<double> y0R,                     // variables values at t=tInit
		         vector<double> parmsR,                  // parameters values
		         vector<double> samplesExogVarR,         // samples for exogenous variables
		         vector<int> nSamplesVarExogVarR,        // number of samples for each exogenous variable
		         int nVarExogVar                         // number of exogenous variables
			 ) {

	
	// Convert vector to ptrs
	double* y0 = &y0R[0];
	double* parms = &parmsR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];
	
	int parmsLength = parmsR.size();
	
	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	//exogVarPeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	modelREuler<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar);
	vector<double> outR(model.getNRowOutSolve());
	double* out = &outR[0];
	model.solve(y0, parms, out);
	return outR;
}

//For now, cubic interpolation is imposed by default (in R)
// [[Rcpp::export]]
vector<double> RK4FixedForR(int nt,                                 // number of points to compute for trajectory
		            double tInit,				 // time at beginning of period to simulate
			    double tEnd,				 // time at end of period to simulate
			    int nV,                                 // number of variables of the system
			    int nIV,                                // number of exogenous variables
		            vector<double> y0R,                     // variables values at t=tInit
		            vector<double> parmsR,                  // parameters values
		            vector<double> samplesExogVarR,         // samples for exogenous variables
		            vector<int> nSamplesVarExogVarR,        // number of samples for each exogenous variable
		            int nVarExogVar                         // number of exogenous variables
			    ) {
	
	// Convert vector to ptrs
	double* y0 = &y0R[0];
	double* parms = &parmsR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];
	
	int parmsLength = parmsR.size();
	
	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	//exogVarPeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	modelRRK4<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar);
	vector<double> outR(model.getNRowOutSolve());
	double* out = &outR[0];
	model.solve(y0, parms, out);
	return outR;
}

//For now, cubic interpolation is imposed by default (in R)
// [[Rcpp::export]]
vector<double> dopriForR(int nt,                                 // number of points to compute for trajectory
		         double tInit,				 // time at beginning of period to simulate
			 double tEnd,				 // time at end of period to simulate
			 int nV,                                 // number of variables of the system
			 int nIV,                                // number of exogenous variables
		         vector<double> y0R,                     // variables values at t=tInit
		         vector<double> parmsR,                  // parameters values
		         vector<double> samplesExogVarR,         // samples for exogenous variables
		         vector<int> nSamplesVarExogVarR,        // number of samples for each exogenous variable
		         int nVarExogVar,                        // number of exogenous variables
			 Rcpp::List argsDopri   		 // List of arguments for variable step method
			 ) {
	double atol = argsDopri["atol"];
	double rtol = argsDopri["rtol"];
	double fac = argsDopri["fac"];
	double facMin = argsDopri["facMin"];
	double facMax = argsDopri["facMax"];
	int nStepMax = argsDopri["nStepMax"];
	double hInit = argsDopri["hInit"];
	double hMin = argsDopri["hMin"];
        double hMax = argsDopri["hMax"];	

	
	// Convert vector to ptrs
	double* y0 = &y0R[0];
	double* parms = &parmsR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];
	
	int parmsLength = parmsR.size();
	
	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	//exogVarPeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	modelRDopri<double> model(nV, nIV, tInit, tEnd, nt, 
			          atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax, 
				  &myExogVar);
	vector<double> outR(model.getNRowOutSolve());
	double* out = &outR[0];
	model.solve(y0, parms, out);
	return outR;
}

//For now, cubic interpolation is imposed by default (in R)
// [[Rcpp::export]]
vector<double> dopriStiffForR(int nt,                                 // number of points to compute for trajectory
							  double tInit,				 // time at beginning of period to simulate
							  double tEnd,				 // time at end of period to simulate
							  int nV,                                 // number of variables of the system
							  int nIV,                                // number of exogenous variables
							  vector<double> y0R,                     // variables values at t=tInit
							  vector<double> parmsR,                  // parameters values
							  vector<double> samplesExogVarR,         // samples for exogenous variables
							  vector<int> nSamplesVarExogVarR,        // number of samples for each exogenous variable
							  int nVarExogVar,                        // number of exogenous variables
						      Rcpp::List argsDopri   		 // List of arguments for variable step method
			 ) {
	double atol = argsDopri["atol"];
	double rtol = argsDopri["rtol"];
	double fac = argsDopri["fac"];
	double facMin = argsDopri["facMin"];
	double facMax = argsDopri["facMax"];
	int nStepMax = argsDopri["nStepMax"];
	double hInit = argsDopri["hInit"];
	double hMin = argsDopri["hMin"];
    double hMax = argsDopri["hMax"];
    int nStiffMax = argsDopri["nStiffMax"];
    int nStiffSuccessiveMax = argsDopri["nStiffSuccessiveMax"];

	// Convert vector to ptrs
	double* y0 = &y0R[0];
	double* parms = &parmsR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];

	int parmsLength = parmsR.size();

	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);

	modelRDopriStiff<double> model(nV, nIV, tInit, tEnd, nt,
			          atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax, nStiffMax, nStiffSuccessiveMax,
				  &myExogVar);
	vector<double> outR(model.getNRowOutSolve());
	double* out = &outR[0];
	model.solve(y0, parms, out);
	return outR;
}

// [[Rcpp::export]]
vector<double> NRForR(int nt,
		      double tInit,
		      double tEnd,
		      int nV, 
		      int nIV,
		      vector<double> y0R,
		      vector<double> parmsR, 
		      vector<double> samplesExogVarR,   
		      vector<int> nSamplesVarExogVarR,
		      int nVarExogVar, 
		      double epsilonJacobian, 
		      double tol, 
		      int nIterMax) {
	
	// Convert vector to ptrs
	double* y0 = &y0R[0];
	double* parms = &parmsR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];
	
	int parmsLength = parmsR.size();	
	
	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	
	// define model and search equilibrium position
	modelREuler<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar); // Only needed to call Func(), Euler's method not used
	NR<double> myNewtonRaphson(&model, epsilonJacobian, tol, nIterMax);
	myNewtonRaphson.NRMethod(tInit, y0, parms);	
	
	return y0R;
}

// [[Rcpp::export]]
vector<double> SAEulerForR(int nt, double tInit, double tEnd,
						   int nV, int nIV, vector<double> y0R,
						   vector<double> parmsDefaultR, vector<double> samplesExogVarR,
						   vector<int> nSamplesVarExogVarR, int nVarExogVar,
						   int nParmsSetSA, int nParmsSA, vector<int> parmsPosSAR,
						   vector<double> allParmsSAR, bool fullTrajectory) {
	//
	// Convert vector to ptrs
	double* y0 = &y0R[0];
	double* parmsDefault = &parmsDefaultR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];
	int parmsLength = parmsDefaultR.size();	
	int* parmsPosSA = &parmsPosSAR[0];
	double* allParmsSA = &allParmsSAR[0];
	
	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	
	// Define model and run simulation
	modelREuler<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar);
	SA<double> mySA(&model, parmsDefault, parmsLength, 
			nParmsSetSA, nParmsSA, parmsPosSA, allParmsSA); 
	if(fullTrajectory==true) {
		vector<double> outR(model.getNRowOutSolve()*nParmsSetSA);
		double* out= &outR[0];
		mySA.sensitivityAnalysisFullTrajectory(y0, out);
		return outR;
	}
	else {
		vector<double> outR(model.getNRowOut()*nParmsSetSA);
		double* out= &outR[0];
		mySA.sensitivityAnalysisLastPoint(y0, out);
		return outR;
	}
}

// [[Rcpp::export]]
vector<double> SARK4FixedForR(int nt, double tInit, double tEnd,
							  int nV, int nIV, vector<double> y0R,
							  vector<double> parmsDefaultR, vector<double> samplesExogVarR,
							  vector<int> nSamplesVarExogVarR, int nVarExogVar,
							  int nParmsSetSA, int nParmsSA, vector<int> parmsPosSAR,
							  vector<double> allParmsSAR, bool fullTrajectory) {
	//
	// Convert vector to ptrs
	double* y0 = &y0R[0];
	double* parmsDefault = &parmsDefaultR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];
	int parmsLength = parmsDefaultR.size();	
	int* parmsPosSA = &parmsPosSAR[0];
	double* allParmsSA = &allParmsSAR[0];
	
	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	
	// Define model and run simulation
	modelRRK4<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar); 
	SA<double> mySA(&model, parmsDefault, parmsLength, 
			nParmsSetSA, nParmsSA, parmsPosSA, allParmsSA); 
	if(fullTrajectory==true) {
		vector<double> outR(model.getNRowOutSolve()*nParmsSetSA);
		double* out= &outR[0];
		mySA.sensitivityAnalysisFullTrajectory(y0, out);
		return outR;
	}
	else {
		vector<double> outR(model.getNRowOut()*nParmsSetSA);
		double* out= &outR[0];
		mySA.sensitivityAnalysisLastPoint(y0, out);
		return outR;
	}
}

// [[Rcpp::export]]
vector<double> SADopriForR(int nt, double tInit, double tEnd,
			   	   	   	   int nV, int nIV, vector<double> y0R,
						   vector<double> parmsDefaultR, vector<double> samplesExogVarR,
						   vector<int> nSamplesVarExogVarR, int nVarExogVar,
						   Rcpp::List argsDopri,
						   int nParmsSetSA, int nParmsSA, vector<int> parmsPosSAR,
						   vector<double> allParmsSAR, bool fullTrajectory) {
	// init parms for Dopri
	double atol = argsDopri["atol"];
	double rtol = argsDopri["rtol"];
	double fac = argsDopri["fac"];
	double facMin = argsDopri["facMin"];
	double facMax = argsDopri["facMax"];
	int nStepMax = argsDopri["nStepMax"];
	double hInit = argsDopri["hInit"];
	double hMin = argsDopri["hMin"];
        double hMax = argsDopri["hMax"];	
	
	// Convert vector to ptrs
	double* y0 = &y0R[0];
	double* parmsDefault = &parmsDefaultR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];
	int parmsLength = parmsDefaultR.size();	
	int* parmsPosSA = &parmsPosSAR[0];
	double* allParmsSA = &allParmsSAR[0];
	
	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	
	// Define model and run simulation
	modelRDopri<double> model(nV, nIV, tInit, tEnd, nt, atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax, &myExogVar);
	SA<double> mySA(&model, parmsDefault, parmsLength, 
			nParmsSetSA, nParmsSA, parmsPosSA, allParmsSA); 
	if(fullTrajectory==true) {
		vector<double> outR(model.getNRowOutSolve()*nParmsSetSA);
		double* out= &outR[0];
		mySA.sensitivityAnalysisFullTrajectory(y0, out);
		return outR;
	}
	else {
		vector<double> outR(model.getNRowOut()*nParmsSetSA);
		double* out= &outR[0];
		mySA.sensitivityAnalysisLastPoint(y0, out);
		return outR;
	}
}

// [[Rcpp::export]]
vector<double> SADopriStiffForR(int nt, double tInit, double tEnd,
								int nV, int nIV, vector<double> y0R,
								vector<double> parmsDefaultR, vector<double> samplesExogVarR,
								vector<int> nSamplesVarExogVarR, int nVarExogVar,
								Rcpp::List argsDopri,
								int nParmsSetSA, int nParmsSA, vector<int> parmsPosSAR,
								vector<double> allParmsSAR, bool fullTrajectory) {
	// init parms for Dopri
	double atol = argsDopri["atol"];
	double rtol = argsDopri["rtol"];
	double fac = argsDopri["fac"];
	double facMin = argsDopri["facMin"];
	double facMax = argsDopri["facMax"];
	int nStepMax = argsDopri["nStepMax"];
	double hInit = argsDopri["hInit"];
	double hMin = argsDopri["hMin"];
    double hMax = argsDopri["hMax"];
    int nStiffMax = argsDopri["nStiffMax"];
    int nStiffSuccessiveMax = argsDopri["nStiffSuccessiveMax"];
	// Convert vector to ptrs
	double* y0 = &y0R[0];
	double* parmsDefault = &parmsDefaultR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];
	int parmsLength = parmsDefaultR.size();
	int* parmsPosSA = &parmsPosSAR[0];
	double* allParmsSA = &allParmsSAR[0];

	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);

	// Define model and run simulation
	modelRDopriStiff<double> model(nV, nIV, tInit, tEnd, nt, atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax, nStiffMax, nStiffSuccessiveMax, &myExogVar);
	SA<double> mySA(&model, parmsDefault, parmsLength,
			nParmsSetSA, nParmsSA, parmsPosSA, allParmsSA);
	if(fullTrajectory==true) {
		vector<double> outR(model.getNRowOutSolve()*nParmsSetSA);
		double* out= &outR[0];
		mySA.sensitivityAnalysisFullTrajectory(y0, out);
		return outR;
	}
	else {
		vector<double> outR(model.getNRowOut()*nParmsSetSA);
		double* out= &outR[0];
		mySA.sensitivityAnalysisLastPoint(y0, out);
		return outR;
	}
}


// [[Rcpp::export]]
vector<int> basinBruteForceForR(int nt,	double tInit, double tEnd,
				int nV, int nIV,
				vector<double> yEqR, vector<double> parmsR, 
				vector<double> samplesExogVarR, vector<int> nSamplesVarExogVarR, int nVarExogVar, 
				bool useDist,
				double tol, std::vector<double> boundsSDR, std::vector<bool> useForConvR,
				std::vector<double> gridR, int gridSize, int solver, Rcpp::List argsDopri) {
	// init parms for Dopri
	double atol = argsDopri["atol"];
	double rtol = argsDopri["rtol"];
	double fac = argsDopri["fac"];
	double facMin = argsDopri["facMin"];
	double facMax = argsDopri["facMax"];
	int nStepMax = argsDopri["nStepMax"];
	double hInit = argsDopri["hInit"];
	double hMin = argsDopri["hMin"];
    double hMax = argsDopri["hMax"];
    int nStiffMax = argsDopri["nStiffMax"];
    int nStiffSuccessiveMax = argsDopri["nStiffSuccessiveMax"];


	// Convert vector to ptrs
	double* parms = &parmsR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];
	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);

	double* yEq = &yEqR[0];
	double* boundsSD = &boundsSDR[0];
	double* grid = &gridR[0];
	// Need deep copy for vector<bool> due to stdlib implementation
	bool useForConv[nV];
	for(int it=0; it<nV; it++) useForConv[it] = useForConvR[it];	

	std::vector<int> outR(gridSize);

	convModel<double> myConv = {useDist, useForConv, yEq, tol, boundsSD};
	gridModel<double> myGrid = {grid, gridSize};

	// define model and run simulation
	if (solver==0) { //dopri
		modelRDopri<double> model(nV, nIV, tInit, tEnd, nt, atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax, &myExogVar);
		basin<double> myBasin(&model, parms, myGrid, myConv);
		myBasin.computeBasin(&outR[0]);
	} else if (solver==1) { //RK4Fixed
		modelRRK4<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar);
		basin<double> myBasin(&model, parms, myGrid, myConv);
		myBasin.computeBasin(&outR[0]);
	} else if (solver==2) { //Euler
		modelREuler<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar);
		basin<double> myBasin(&model, parms, myGrid, myConv);
		myBasin.computeBasin(&outR[0]);
	} else if (solver==3) {
		modelRDopriStiff<double> model(nV, nIV, tInit, tEnd, nt, atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax, nStiffMax, nStiffSuccessiveMax, &myExogVar);
		basin<double> myBasin(&model, parms, myGrid, myConv);
		myBasin.computeBasin(&outR[0]);
	}
	return outR;
}

// [[Rcpp::export]]
Rcpp::List CMAESForR(int nt,
		     double tInit,
		     double tEnd,
		     int nV,
		     int nIV,
		     vector<double> parmsR, 
		     vector<double> samplesExogVarR,   
		     vector<int> nSamplesVarExogVarR,
		     int nVarExogVar, 
		     int nVarMinDist,
		     vector<int> varMinDistR,
		     vector<double> dataMinDistR,
		     vector<double> samplingTimeMinDistR,   
		     vector<double> pointsWeightMinDistR,
		     vector<int>  nObsVarMinDistR, 
		     vector<bool> parmsToOptimizeR,
		     vector<double> parmsLowerR, 
		     vector<double> parmsUpperR,
		     int lambda, 
		     double sigma, 
		     int nIterMax, 
		     double tolCMAES, 
		     bool normalizeParms, 
			 int nReDrawMax,
		     int solver, 
		     Rcpp::List argsDopri) {
	// init parms for Dopri
	double atol = argsDopri["atol"];
	double rtol = argsDopri["rtol"];
	double fac = argsDopri["fac"];
	double facMin = argsDopri["facMin"];
	double facMax = argsDopri["facMax"];
	int nStepMax = argsDopri["nStepMax"];
	double hInit = argsDopri["hInit"];
	double hMin = argsDopri["hMin"];
    double hMax = argsDopri["hMax"];
    int nStiffMax = argsDopri["nStiffMax"];
    int nStiffSuccessiveMax = argsDopri["nStiffSuccessiveMax"];

	int parmsLength = parmsR.size();	
	
	// Convert vector to ptrs
	double* parms = &parmsR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];
	int* varMinDist = &varMinDistR[0];
	double* dataMinDist = &dataMinDistR[0];
	double* samplingTimeMinDist = &samplingTimeMinDistR[0];
	double* pointsWeightMinDist = &pointsWeightMinDistR[0];
	int* nObsVarMinDist = &nObsVarMinDistR[0];
	double* parmsLower = &parmsLowerR[0];
	double* parmsUpper = &parmsUpperR[0];
	bool parmsToOptimize[parmsToOptimizeR.size()];
	int nParmsToOptimize=0;
	for (int it=0; it<parmsLength; it++) {
		parmsToOptimize[it] = parmsToOptimizeR[it];
		if (parmsToOptimize[it]==true) nParmsToOptimize++;
	}
	
	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	//exogVarPeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	
	// Define vector of non-fixed vectors for optimisation
	int indexOut=0;
	arma::vec outR(nParmsToOptimize);
	for (int it=0; it<parmsLength; it++) {
		if(parmsToOptimize[it]==true) {
			outR(indexOut) = parms[it];
			indexOut++;
		}
	}

	double dist = 0;

	if (solver==0) {
		modelRDopri<double> model(nV, nIV, tInit, tEnd, nt, atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax, &myExogVar);
		minDist<double> myMinDist(&model, parmsLength, nV,
			  	    	  nVarMinDist, varMinDist, dataMinDist, samplingTimeMinDist, pointsWeightMinDist, nObsVarMinDist);
		fixedObjFunc<double, minDist<double>> myFixedObjFunc(myMinDist, parms, parmsToOptimize);
		CMAES<double, fixedObjFunc<double, minDist<double>>> myCMAES(myFixedObjFunc,
		       				    	    lambda, sigma, nIterMax, tolCMAES,
						    	    parmsLower, parmsUpper, nReDrawMax, UseParallel);
		dist = myCMAES.Optimize(outR);

	} else if (solver==1) {
		modelRRK4<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar);
		minDist<double> myMinDist(&model, parmsLength, nV,
			  	    	  nVarMinDist, varMinDist, dataMinDist, samplingTimeMinDist, pointsWeightMinDist, nObsVarMinDist);
		fixedObjFunc<double, minDist<double>> myFixedObjFunc(myMinDist, parms, parmsToOptimize);
		CMAES<double, fixedObjFunc<double, minDist<double>>> myCMAES(myFixedObjFunc,
		       				    	    lambda, sigma, nIterMax, tolCMAES,
						    	    parmsLower, parmsUpper, nReDrawMax, UseParallel);
		dist = myCMAES.Optimize(outR);

	} else if (solver==2) {
		modelREuler<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar);
		minDist<double> myMinDist(&model, parmsLength, nV,
			  	    	  nVarMinDist, varMinDist, dataMinDist, samplingTimeMinDist, pointsWeightMinDist, nObsVarMinDist);
		fixedObjFunc<double, minDist<double>> myFixedObjFunc(myMinDist, parms, parmsToOptimize);
		CMAES<double, fixedObjFunc<double, minDist<double>>> myCMAES(myFixedObjFunc,
		       				    	    							 lambda, sigma, nIterMax, tolCMAES,
																	 parmsLower, parmsUpper, nReDrawMax, UseParallel);
		dist = myCMAES.Optimize(outR);
	} else if (solver==3) {
		modelRDopriStiff<double> model(nV, nIV, tInit, tEnd, nt, atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax, nStiffMax, nStiffSuccessiveMax, &myExogVar);
		minDist<double> myMinDist(&model, parmsLength, nV,
			  	    	  nVarMinDist, varMinDist, dataMinDist, samplingTimeMinDist, pointsWeightMinDist, nObsVarMinDist);
		fixedObjFunc<double, minDist<double>> myFixedObjFunc(myMinDist, parms, parmsToOptimize);
		CMAES<double, fixedObjFunc<double, minDist<double>>> myCMAES(myFixedObjFunc,
																	 lambda, sigma, nIterMax, tolCMAES,
																	 parmsLower, parmsUpper, nReDrawMax, UseParallel);
		dist = myCMAES.Optimize(outR);
	}
	return Rcpp::List::create(Rcpp::Named("parms") = outR, 
				  Rcpp::Named("dist") = dist);
}

/*
// [[Rcpp::export]]
Rcpp::List CMAESLLKForR(int nt,
						double tInit,
						double tEnd,
						int nV,
						int nIV,
						vector<double> y0R,
						vector<double> parmsR,
						vector<double> samplesExogVarR,   
						vector<int> nSamplesVarExogVarR,
						int nVarExogVar, 
						int nVarMinDist,
						vector<int> varMinDistR,
						vector<double> dataMinDistR,
						vector<double> samplingTimeMinDistR,   
						vector<double> pointsWeightMinDistR,
						vector<int>  nObsVarMinDistR, 
						vector<bool> parmsToOptimizeR,
						vector<double> parmsLowerR, 
						vector<double> parmsUpperR,
						int lambda, 
						double sigma, 
						int nIterMax, 
						double tolCMAES, 
						bool normalizeParms, 
						int nReDrawMax,
						int solver, 
						Rcpp::List argsDopri) {
	// init parms for Dopri
	double atol = argsDopri["atol"];
	double rtol = argsDopri["rtol"];
	double fac = argsDopri["fac"];
	double facMin = argsDopri["facMin"];
	double facMax = argsDopri["facMax"];
	int nStepMax = argsDopri["nStepMax"];
	double hInit = argsDopri["hInit"];
	double hMin = argsDopri["hMin"];
    double hMax = argsDopri["hMax"];
    int nStiffMax = argsDopri["nStiffMax"];
    int nStiffSuccessiveMax = argsDopri["nStiffSuccessiveMax"];

	int parmsLength = parmsR.size();	
	
	// Convert vector to ptrs
	double* y0 = &y0R[0];
	double* parms = &parmsR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];
	int* varMinDist = &varMinDistR[0];
	double* dataMinDist = &dataMinDistR[0];
	double* samplingTimeMinDist = &samplingTimeMinDistR[0];
	double* pointsWeightMinDist = &pointsWeightMinDistR[0];
	int* nObsVarMinDist = &nObsVarMinDistR[0];
	double* parmsLower = &parmsLowerR[0];
	double* parmsUpper = &parmsUpperR[0];
	bool parmsToOptimize[parmsToOptimizeR.size()];
	int nParmsToOptimize=0;
	for (int it=0; it<parmsLength; it++) {
		parmsToOptimize[it] = parmsToOptimizeR[it];
		if (parmsToOptimize[it]==true) nParmsToOptimize++;
	}
	
	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	//exogVarPeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	
	// Define vector of non-fixed vectors for optimisation
	int indexOut=0;
	arma::vec outR(nParmsToOptimize);
	for (int it=0; it<parmsLength; it++) {
		if(parmsToOptimize[it]==true) {
			outR(indexOut) = parms[it];
			indexOut++;
		}
	}

	double dist = 0;

	if (solver==0) {
		modelRDopri<double> model(nV, nIV, tInit, tEnd, nt, atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax, &myExogVar);
		minDist<double> myMinDist(&model, parmsLength, y0,
			  	    	  nVarMinDist, varMinDist, dataMinDist, samplingTimeMinDist, pointsWeightMinDist, nObsVarMinDist);
		fixedObjFunc<double, minDist<double>> myFixedObjFunc(myMinDist, parms, parmsToOptimize);
		CMAES<double, fixedObjFunc<double, minDist<double>>> myCMAES(myFixedObjFunc,
		       				    	    lambda, sigma, nIterMax, tolCMAES,
						    	    parmsLower, parmsUpper, nReDrawMax, UseParallel);
		dist = myCMAES.Optimize(outR);

	} else if (solver==1) {
		modelRRK4<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar);
		minDist<double> myMinDist(&model, parmsLength, y0,
			  	    	  nVarMinDist, varMinDist, dataMinDist, samplingTimeMinDist, pointsWeightMinDist, nObsVarMinDist);
		fixedObjFunc<double, minDist<double>> myFixedObjFunc(myMinDist, parms, parmsToOptimize);
		CMAES<double, fixedObjFunc<double, minDist<double>>> myCMAES(myFixedObjFunc,
		       				    	    lambda, sigma, nIterMax, tolCMAES,
						    	    parmsLower, parmsUpper, nReDrawMax, UseParallel);
		dist = myCMAES.Optimize(outR);

	} else if (solver==2) {
		modelREuler<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar);
		minDist<double> myMinDist(&model, parmsLength, y0,
			  	    	  nVarMinDist, varMinDist, dataMinDist, samplingTimeMinDist, pointsWeightMinDist, nObsVarMinDist);
		fixedObjFunc<double, minDist<double>> myFixedObjFunc(myMinDist, parms, parmsToOptimize);
		CMAES<double, fixedObjFunc<double, minDist<double>>> myCMAES(myFixedObjFunc,
		       				    	    							 lambda, sigma, nIterMax, tolCMAES,
																	 parmsLower, parmsUpper, nReDrawMax, UseParallel);
		dist = myCMAES.Optimize(outR);
	} else if (solver==3) {
		modelRDopriStiff<double> model(nV, nIV, tInit, tEnd, nt, atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax, nStiffMax, nStiffSuccessiveMax, &myExogVar);
		minDist<double> myMinDist(&model, parmsLength, y0,
			  	    	  nVarMinDist, varMinDist, dataMinDist, samplingTimeMinDist, pointsWeightMinDist, nObsVarMinDist);
		fixedObjFunc<double, minDist<double>> myFixedObjFunc(myMinDist, parms, parmsToOptimize);
		CMAES<double, fixedObjFunc<double, minDist<double>>> myCMAES(myFixedObjFunc,
																	 lambda, sigma, nIterMax, tolCMAES,
																	 parmsLower, parmsUpper, nReDrawMax, UseParallel);
		dist = myCMAES.Optimize(outR);
	}
	return Rcpp::List::create(Rcpp::Named("parms") = outR, 
				  Rcpp::Named("dist") = dist);
}
*/



// [[Rcpp::export]]
arma::mat SCMNR(int nt,
	        double tInit,
		double tEnd,
		int nV,
		int nIV,
	   	vector<double> parmsR,
	   	vector<double> samplesExogVarR,
	   	vector<int> nSamplesVarExogVarR,
	   	int nVarExogVar,
	   	vector<int> nCellsR,
	   	vector<double> yLowerR,
	   	vector<double> yUpperR,
	   	vector<double> byCellsR,
	   	double epsilon,
	   	double tolNR,
	   	double nIterMax,
	  	bool postProcessing) {

	// Convert vector to ptrs
	double* parms = &parmsR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];
	
	int parmsLength = parmsR.size();	
	
	int* nCells = &nCellsR[0];
	double* yLower = &yLowerR[0];
	double* yUpper = &yUpperR[0];
	double* byCells = &byCellsR[0];

	// DEFINE CLASSES FOR COMPUTATION
	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	modelREuler<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar);
	cellMapNR<double> myCellMapNR(&model, parms, nCells, byCells, yLower, yUpper, epsilon, tolNR, nIterMax);
	SCM<cellMapNR<double>> mySCMNR(myCellMapNR, false);
	
	// SIMULATION
	mySCMNR.computeSCM(postProcessing);

	// OUTPUT FORMATTING FOR R AND POSTPROCESSING	
	std::vector<std::vector<int>> allAttractors = *mySCMNR.returnAttractors();
	arma::mat out(allAttractors.size(), nV);
	for (size_t it=0; it<allAttractors.size(); it++) {
		double tempY[nV];
	        mySCMNR.cellToCoord(allAttractors[it][0], tempY);
		if (postProcessing==true) myCellMapNR.NRMethod(tInit, tempY, parms);
		for (size_t j=0; j<nV; j++) {
			out(it, j) = tempY[j];
		}
	}
	return out;

}

// [[Rcpp::export]]
arma::mat HCMSubdivisionNR(int nt,
			   double tInit,
			   double tEnd, 
			   int nV, 
			   int nIV,
			   vector<double> parmsR, 
			   vector<double> samplesExogVarR,   
			   vector<int> nSamplesVarExogVarR,
			   int nVarExogVar, 
			   vector<int> nCellsR, 
			   vector<double> yLowerR, 
			   vector<double> yUpperR, 
			   vector<double> byCellsR, 
			   double epsilon, 
			   double tolNR, 
			   double nIterMax, 
			   int nSubdivisionMax, 
			   int nSubdivisionGCM, 
			   std::vector<double> gridCellR, 
			   int nEltsGridCell, 
			   bool postProcessing) {
	
	// Convert vector to ptrs
	double* parms = &parmsR[0];
	double* samplesExogVar = &samplesExogVarR[0];
	int* nSamplesVarExogVar = &nSamplesVarExogVarR[0];
	double* gridCell = &gridCellR[0];
	
	int parmsLength = parmsR.size();
	
	int* nCells = &nCellsR[0];
	double* yLower = &yLowerR[0];
	double* yUpper = &yUpperR[0];
	double* byCells = &byCellsR[0];
	
	// DEFINE CLASSES FOR COMPUTATION
	exogVarCubicSplinePeriodic2<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	modelREuler<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar);
	cellMapSubdivisionNR<double, TInt> myCellMapSubNR(&model, parms, 
				   		          nCells, byCells, 
						          yLower, yUpper, epsilon, tolNR, nIterMax, 
						          gridCell, nEltsGridCell);
	HCMSubdivision<cellMapSubdivisionNR<double, TInt>, TInt> myHCMSubNR(myCellMapSubNR, nSubdivisionMax, nSubdivisionGCM);

	// SIMULATION
	myHCMSubNR.computeHCMSubdivision();

	std::vector<TInt> allAttractors = *myHCMSubNR.returnResults();
	arma::mat out(allAttractors.size(), nV);
	for (size_t i=0; i<allAttractors.size(); i++) {
		double tempY[nV];
	        myHCMSubNR.cellToCoord(allAttractors[i], tempY);
		if (postProcessing==true) myCellMapSubNR.NRMethod(tInit, tempY, parms);
		for (size_t j=0; j<nV; j++) {
			out(i, j) = tempY[j];
		}
	}
	return out;
}


