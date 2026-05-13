#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "src/preproc.h"
#include "src/goodwinKeenEuler.h"

using namespace std;


int main() {
	// INIT MAIN PARAMETERS
	int nV = 3;
	int nIV = 5;
	const int parmsLength = 10;
	double parms[parmsLength] = {0.025, 0.02, 0.01, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};	
	int nt = 30001;
	double tInit = 0;
	double tEnd = 300;
	/**************************/
	/**************************/
	/* RUN SIMULATION VANILLA */
	/**************************/
	/**************************/
	cout<<"EULER GK VANILLA"<<endl;	
	// INIT MODEL+SOLVER
	double y0[3] = {0.8364, 0.9674, 0.07109};	
	GoodwinKeenEulerVanilla<double> GKEV(nV, nIV, tInit, tEnd, nt);
	double outLastPoint[GKEV.getNRowOut()];
	double outFull[GKEV.getNRowOutSolve()];
	GKEV.solveLastPoint(y0, parms, outLastPoint);
	GKEV.solve(y0, parms, outFull);

/*        for (int it=0; it<nt; it++) {
		cout<<outFull[it*GKEV.nRowOut()]<<endl;
	} */
	cout<<"init"<<y0[0]<<" "<<y0[1]<<" "<<y0[2]<<endl;
	cout<<"last point: "<<outLastPoint[0]
		<<" "<<outLastPoint[1]
		<<" "<<outLastPoint[2]<<endl;


	/*******************************/
	/*******************************/
	/* RUN SIMULATION WITH EXOGVAR */
	/*******************************/
	/*******************************/
	cout<<"EULER GK EXOG VARIABLE"<<endl;	
	// INIT MODEL+SOLVER
	// Define exogVar
	int nVarExogVar = 2;
	int nSamplesVarExogVar[2] = {11, 3}; 
	double samplesExogVar[14] = {0,1,2,3,4,5,6,7,8,9,10,0, 50, 100};
	
	exogVarPeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	GoodwinKeenEulerEV<double> GKEEV(nV, nIV, tInit, tEnd, nt, &myExogVar);
	double outLastPointEV[GKEEV.getNRowOut()];
	double outFullEV[GKEEV.getNRowOutSolve()];

	GKEEV.solveLastPoint(y0, parms, outLastPointEV);

	//GKEEV.solve(y0, parms, outFullEV);
        for (int it=0; it<nt; it++) {
	//	cout<<outFullEV[it*GKEEV.nRowOut()]<<endl;
	}
	cout<<"init"<<y0[0]<<" "<<y0[1]<<" "<<y0[2]<<endl;
	cout<<"last point: "<<outLastPointEV[0]
		<<" "<<outLastPointEV[1]
		<<" "<<outLastPointEV[2]<<endl;

	cout<<"last point Exog Var: "<<outLastPointEV[6]<<" "<<outLastPointEV[7]<<endl;
	return 0;
}
