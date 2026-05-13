#define UseRCPP 1 		// true if the code is called from R(via RCPP), false if it is directly called from c++


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

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//#include "armadillo"
#include "src/preprocRCPP_R.h"
//#include "src/preproc.h"

#include "src/exogenousVariables.h"
#include "src/modelR.h"
#include "src/solverModelR.h"
#include "src/ODE.h"
#include "src/euler.h"
#include "src/fixedMinDist.h"
#include "src/minDist.h"
#include "src/CMAES.h"

// [[Rcpp::plugins("cpp17")]]

// [[Rcpp::export]]
int temp() {
	// INIT MAIN PARAMETERS
	const int parmsLength = 10;
	double parms[parmsLength] = {0.025, 0.02, 0.01, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};	
	const double tInit = 0.0;
	const double tEnd = 300.0;
	const int nt = 30001;
	
	// Define exogVar
	int nVarExogVar = 2;
	int nSamplesVarExogVar[2] = {11, 3}; 
	double samplesExogVar[14] = {0,1,2,3,4,5,6,7,8,9,10,0, 50, 100};

	const int nV = 3;
	const int nIV = 5;


	/******************************/	
	/******************************/	
	/* CODE DISTANCE MINIMIZATION */
	/******************************/	
	/******************************/	

	const int nVarMinDist = 3;
	int varMinDist[3] = {0, 1, 2};
	double dataMinDist[3*301];
	double samplingTimeMinDist[3*301];
	double pointsWeightMinDist[3*301];
	int nObsVarMinDist[nVarMinDist] = {301, 301, 301};
	for (int it=0; it<301; it++) {
		dataMinDist[it] = 0.83605402;
	        dataMinDist[301+it] = 0.96861321;
	        dataMinDist[602+it] = 0.07019744;
		samplingTimeMinDist[it] = it;
		samplingTimeMinDist[301+it] = it;
		samplingTimeMinDist[602+it] = it;
		pointsWeightMinDist[it] = 1;
		pointsWeightMinDist[301+it] = 1;
		pointsWeightMinDist[602+it] = 1;
	}
	
	double y0[nV] = {0.83605402, 0.96861321, 0.07019744};


	bool parmsToOptimize[parmsLength] = {true, true, true, false, false, false, false, false, false, false};
	double parmsLower[parmsLength] = {0.010, 0.01, 0.0, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};			
	double parmsUpper[parmsLength] = {0.035, 0.03, 0.02, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};		
	const int lambda = 25;
	const double sigma = 0.1;
	const int nIterMax = 200;
	const double tolCMAES = 1e-4;
		
	
	exogVarCubicSplinePeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	modelREuler<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar);
	minDist<double> myMinDist(&model, parmsLength, y0,
		  	    	  nVarMinDist, varMinDist, dataMinDist, samplingTimeMinDist, pointsWeightMinDist, nObsVarMinDist);
	fixedMinDist<double> myFixedMinDist(myMinDist, parms, parmsToOptimize);	      
	CMAES<double, fixedMinDist<double>> myCMAES(myFixedMinDist,
	       				    	    lambda, sigma, nIterMax, tolCMAES, 
					    	    parmsLower, parmsUpper);
	arma::vec parmsPartForTest1(3);
	parmsPartForTest1 = {0.011, 0.012, 0.012};
	double dist = myCMAES.Optimize(parmsPartForTest1);
	cout<<parmsPartForTest1[0]<<" "<<parmsPartForTest1[1]<<" "<<parmsPartForTest1[2]<<endl;
	
	return 0;
}


int main() {
	temp();
	return 0;
}
