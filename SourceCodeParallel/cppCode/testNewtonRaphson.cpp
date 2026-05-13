#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

#include "src/preproc.h"
#include "src/newtonRaphson.h"
#include "src/goodwinKeenRK4.h"
using namespace std;


int main() {
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

	exogVarPeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	const int nV = 3;
	const int nIV = 5;
	GoodwinKeenRK4EV<double> keenModel(nV, nIV, tInit, tEnd, nt, &myExogVar); 
	
	// Parms for Newton-Raphson
	double epsilon = 1e-6;
	double tol = 1e-6;
	int nIterMax = 100;

	NR<double> myNR(&keenModel, epsilon, tol, nIterMax);

	double y0[nV] = {0.7, 0.95, 0.05};

	myNR.NRMethod(tInit, y0, parms);		
	cout<<"equilibrium: "<<y0[0]<<" "<<y0[1]<<" "<<y0[2]<<endl;
	return 0;
}
