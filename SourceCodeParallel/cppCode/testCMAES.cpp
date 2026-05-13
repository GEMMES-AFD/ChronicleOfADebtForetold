#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <cblas.h>

#define ARMA_WARN_LEVEL 1
#include <armadillo>

#include "src/goodwinKeenRK4.h"
#include "src/preproc.h"
#include "src/CMAES.h"
#include "src/fixedObjFunc.h"
#include "src/minDist.h"
using namespace std;


int main() {
	// INIT MAIN PARAMETERS
	const int parmsLength = 10;
	double parms[parmsLength] = {0.025, 0.02, 0.01, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};	
	const double tInit = 0.0;
	const double tEnd = 300.0;
	const int nt = 30001;
	int nReDrawMax = 100;
	// Define exogVar
	int nVarExogVar = 2;
	int nSamplesVarExogVar[2] = {11, 3}; 
	double samplesExogVar[14] = {0,1,2,3,4,5,6,7,8,9,10,0, 50, 100};

	exogVarPeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	const int nV = 3;
	const int nIV = 5;

	GoodwinKeenRK4EV<double> model(nV, nIV, tInit, tEnd, nt, &myExogVar); 

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
	
	double yInit[nV] = {0.83605402, 0.96861321, 0.07019744};


	minDist<double> myMinDist(&model, parmsLength, yInit,
				    nVarMinDist, varMinDist, dataMinDist, samplingTimeMinDist, pointsWeightMinDist, nObsVarMinDist);
	bool parmsToOptimize[parmsLength] = {true, true, true, false, false, false, false, false, false, false};
	double parmsLower[parmsLength] = {0.015, 0.01, 0.0, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};
	double parmsUpper[parmsLength] = {0.035, 0.03, 0.02, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};
	cout<<"EvaluateMinDist: "<<myMinDist.Evaluate(parms)<<endl;
	fixedMinDist<double> myFixedMinDist(myMinDist, parms, parmsToOptimize);	      
	
	const int lambda = 25;
	const double sigma = 0.1;
	const int nIterMax = 200;
	const double tolCMAES = 0.001;
		
	
	// CMAES Class declaration
	CMAES<double, fixedMinDist<double>> myCMAES(myFixedMinDist,
		       				    lambda, sigma, nIterMax, tolCMAES, 
						    parmsLower, parmsUpper, nReDrawMax);

	// Simulation and test	
	arma::vec parmsPartForTest1(3);
	arma::vec parmsPartForTest2(3);
	parmsPartForTest1 = {0.024, 0.019, 0.012};
	parmsPartForTest2 = {0.022, 0.019, 0.012};
	cout<<"Print init parms\n"<<parmsPartForTest1<<endl;
	cout<<parms[0]<<" "<<parms[9]<<endl;
	// To use to test distance estimation, iif parmsPartForTest initialized with values from parmsDouble
	double toPrint1 = myFixedMinDist.Evaluate(parmsPartForTest1);
        //if (standardizeParms==true) myCMAES.standardizeParmsPart(parmsPartForTest);
	double toPrint2 = myCMAES.Evaluate(parmsPartForTest1);
	//if (standardizeParms==true) myCMAES.unstandardizeParmsPart(parmsPartForTest);
	double toPrint3 = myCMAES.Evaluate(parmsPartForTest1);
	double toPrint3bis = myCMAES.Evaluate(parmsPartForTest2);
	//cout<< "test minDist and fixedMinDist give same initial distance: "<<toPrint1<<" "<<toPrint2<<" "<<toPrint3<<endl;
	//cout<<"test EvaluateBounded on two different initial positions: "<<toPrint3<<" "<<toPrint3bis<<endl;
	double toPrint4 = myCMAES.Optimize(parmsPartForTest1);
	cout<<"test(CMAES)\n, dist = "<<toPrint4<<" parms =\n "<<parmsPartForTest1<<endl;	
	









	return 0;
}	
