#include <chrono>

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "src/preproc.h"
#include "src/goodwinKeenDopriStiff.h"

using namespace std;


int main() {
	// INIT MAIN PARAMETERS
	int nV = 3;
	int nIV = 5;
	const int parmsLength = 10;
	double parms[parmsLength] = {0.025, 0.02, 0.01, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};
	int nt = 500;
	double tInit = 0;
	double tEnd = 10;
	cout<<(double)nt/tEnd<<endl;;
//	double y0[3] = {0.8364, 0.9674, 0.07109};

	double y0[3] = {0.85, 0.96, 5};

	// Init parameters for dopri
	double atol   = 1e-4;
	double rtol   = 0;
	double fac    = 0.85;
	double facMin = 0.1;
	double facMax = 4;
	int nStepMax  = 1000;
	double hInit  = 0.01;
	double hMin   = 0.00001;
	double hMax   = 0.1;

	int nStiffMax = 100;
	int nStiffSuccessiveMax = 5;

	/**************************/
	/**************************/
	/* RUN SIMULATION VANILLA */
	/**************************/
	/**************************/
	cout<<"DOPRI GK VANILLA"<<endl;
	// INIT MODEL+SOLVER
	GoodwinKeenDopriStiff<double> GKDopriV(nV, nIV, tInit, tEnd, nt, nullptr, atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax, nStiffMax, nStiffSuccessiveMax);
	double outLastPoint[GKDopriV.getNRowOut()];
	double outFull[GKDopriV.getNRowOutSolve()];

	GKDopriV.solveLastPoint(y0, parms, outLastPoint);
	GKDopriV.solve(y0, parms, outFull);

/*        for (int it=0; it<nt; it++) {
		cout<<outFull[it*GKRK4V.getNRowOut()]<<endl;
	}*/
	cout<<"init"<<y0[0]<<" "<<y0[1]<<" "<<y0[2]<<endl;
	cout<<"last point: "<<outLastPoint[0]
		<<" "<<outLastPoint[1]
		<<" "<<outLastPoint[2]<<endl;
	cout<<"last point: "<<outFull[GKDopriV.getNRowOutSolve()-2*nV-nIV]<<" "<<outFull[GKDopriV.getNRowOutSolve() - 2*nV - nIV+1]<<" "<<outFull[GKDopriV.getNRowOutSolve() - 2*nV - nIV+2]<<endl;

	return 0;
}
