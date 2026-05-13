#include <chrono>	

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "src/preproc.h"
#include "src/goodwinKeenDopri.h"

using namespace std;


int main() {
	// INIT MAIN PARAMETERS
	int nV = 3;
	int nIV = 5;
	const int parmsLength = 10;
	double parms[parmsLength] = {0.025, 0.02, 0.01, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};	
	int nt = 500;
	double tInit = 0;
	double tEnd = 1000;
	cout<<(double)nt/tEnd<<endl;;
//	double y0[3] = {0.8364, 0.9674, 0.07109};	

	double y0[3] = {0.86, 0.98, 0.10};	

	// Init parameters for dopri
	double atol = 1e-9;
	double rtol = 1e-9;
	double fac = 0.85;
	double facMin = 0.1;
	double facMax = 4;
	int nStepMax = 1000;
	double hInit = 0.01;
	double hMin = 0.0001;
        double hMax = 0.2;	

	/**************************/
	/**************************/
	/* RUN SIMULATION VANILLA */
	/**************************/
	/**************************/
	cout<<"DOPRI GK VANILLA"<<endl;	
	// INIT MODEL+SOLVER
	GoodwinKeenDopriVanilla<double> GKDopriV(nV, nIV, tInit, tEnd, nt, nullptr, atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax);
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

	/*******************************/
	/*******************************/
	/* RUN SIMULATION WITH EXOGVAR */
	/*******************************/
	/*******************************/
	cout<<"DOPRI GK EXOG VARIABLE"<<endl;	
	// INIT MODEL+SOLVER
	// Define exogVar
	int nVarExogVar = 2;
	int nSamplesVarExogVar[2] = {11, 3}; 
	double samplesExogVar[14] = {0,1,2,3,4,5,6,7,8,9,10,0, 50, 100};
	double samplingTimeExogVar[14] = {0,1,2,3,4,5,6,7,8,9,10,0, 50, 100};

	exogVarPeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	GoodwinKeenDopriEV<double> GKDopriEV(nV, nIV, tInit, tEnd, nt, &myExogVar, atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax);
	double outLastPointEV[GKDopriEV.getNRowOut()];
	double outFullEV[GKDopriEV.getNRowOutSolve()];

	GKDopriEV.solveLastPoint(y0, parms, outLastPointEV);

	//GKEEV.solve(y0, parms, outFullEV);
        for (int it=0; it<nt; it++) {
	//	cout<<outFullEV[it*GKEEV.getNRowOut()]<<endl;
	}
	cout<<"init"<<y0[0]<<" "<<y0[1]<<" "<<y0[2]<<endl;
	cout<<"last point: "<<outLastPointEV[0]
		<<" "<<outLastPointEV[1]
		<<" "<<outLastPointEV[2]<<endl;

	cout<<"last point Exog Var: "<<outLastPointEV[6]<<" "<<outLastPointEV[7]<<endl;
	cout<<"DONE !"<<endl;

	//measure computation time: 

	int nTrials = 1000;
	
	auto start = std::chrono::steady_clock::now();	
	for (int it=0; it<nTrials; it++) {
		GKDopriEV.solveLastPoint(y0, parms, outLastPointEV);
	}
	auto diff = std::chrono::steady_clock::now() - start;
    	auto outTime = std::chrono::duration_cast<std::chrono::milliseconds>(diff);
	cout<<"computation time: "<<outTime.count()<<endl;

	return 0;
}
