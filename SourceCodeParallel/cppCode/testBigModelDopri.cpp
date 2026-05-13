#include <chrono>	

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "src/preproc.h"
#include "src/bigModelDopri.h"

using namespace std;


int main() {
	// INIT MAIN PARAMETERS
	int nV = 100;
	int nIV = 5;
	const int parmsLength = 10;
	double parms[parmsLength] = {0.025, 0.02, 0.01, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};	
	int nt = 1000;
	double tInit = 0;
	double tEnd = 1000;
	cout<<(double)nt/tEnd<<endl;;
//	double y0[3] = {0.8364, 0.9674, 0.07109};	

	double y0[nV];
	for (int it=0; it<nV; it++) y0[it] = it;	
	
	// Init parameters for dopri
	double atol = 1e-4;
	double rtol = 1e-4;
	double fac = 0.85;
	double facMin = 0.1;
	double facMax = 4;
	int nStepMax = 1000;
	double hInit = 0.01;
	double hMin = 0.0001;
       	double hMax = 0.5;	

	/**************************/
	/**************************/
	/*     RUN SIMULATION     */
	/**************************/
	/**************************/
	// INIT MODEL+SOLVER
	bigModelDopri<double> bigModelDopri(nV, nIV, tInit, tEnd, nt, nullptr, atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax);
	double outLastPoint[bigModelDopri.getNRowOut()];
	double outFull[bigModelDopri.getNRowOutSolve()];

	bigModelDopri.solveLastPoint(y0, parms, outLastPoint);
	bigModelDopri.solve(y0, parms, outFull);

/*        for (int it=0; it<nt; it++) {
		cout<<outFull[it*GKRK4V.getNRowOut()]<<endl;
	}*/
	cout<<"init"<<y0[0]<<" "<<y0[1]<<" "<<y0[2]<<endl;
	cout<<"last point: "<<outLastPoint[0]
		<<" "<<outLastPoint[1]
		<<" "<<outLastPoint[2]<<endl;
	cout<<"last point: "<<outFull[bigModelDopri.getNRowOutSolve()-2*nV-nIV]<<" "<<outFull[bigModelDopri.getNRowOutSolve() - 2*nV - nIV+1]<<" "<<outFull[bigModelDopri.getNRowOutSolve() - 2*nV - nIV+2]<<endl;

	// same but with RK4 for comparison
	double y0RK4[nV];
	for (int it=0; it<nV; it++) y0RK4[it] = it;	
	int ntRK4 = 10000;
	double tInitRK4 = 0;
	double tEndRK4 = 1000;
	bigModelRK4<double> bigModelRK4(nV, nIV, tInitRK4, tEndRK4, ntRK4, nullptr);
	cout<<bigModelRK4.getNRowOut()<<" "<<bigModelRK4.getNRowOutSolve()<<endl;
	double outLastPointRK4[bigModelRK4.getNRowOut()];
	//double outFullRK4[bigModelRK4.getNRowOutSolve()];
	bigModelRK4.solveLastPoint(y0RK4, parms, outLastPointRK4);
	cout<<"init"<<y0[0]<<" "<<y0[1]<<" "<<y0[2]<<endl;
	cout<<"last point: "<<outLastPointRK4[0]
		<<" "<<outLastPointRK4[1]
		<<" "<<outLastPointRK4[2]<<endl;

	//measure computation time: 

	int nTrials = 100;
	
	auto startRK4 = std::chrono::steady_clock::now();	
	for (int it=0; it<nTrials; it++) {
		bigModelRK4.solveLastPoint(y0, parms, outLastPoint);
	}
	auto diffRK4 = std::chrono::steady_clock::now() - startRK4;
    	auto outTimeRK4 = std::chrono::duration_cast<std::chrono::milliseconds>(diffRK4);

	auto startDopri = std::chrono::steady_clock::now();	
	for (int it=0; it<nTrials; it++) {
		bigModelDopri.solveLastPoint(y0, parms, outLastPoint);
	}
	auto diffDopri = std::chrono::steady_clock::now() - startDopri;
    	auto outTimeDopri = std::chrono::duration_cast<std::chrono::milliseconds>(diffDopri);
	cout<<"coucou"<<endl;
	cout<<"computation time with Dopri: "<<outTimeDopri.count()<<endl;
	cout<<"computation time with RK4: "<<outTimeRK4.count()<<endl;

	return 0;
}
