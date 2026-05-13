#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <cblas.h>
#define ARMA_WARN_LEVEL 1
#define ARMA_DONT_USE_WRAPPER 1
#include <armadillo>

#include <chrono>
#include "src/preproc.h"
#include "src/goodwinKeenRK4.h"
#include "src/SCM.h"
#include "src/cellMapNR.h"

using namespace std;


int main() {
	// INIT MAIN PARAMETERS
	const int parmsLength = 10;
	double parms[parmsLength] = {0.025, 0.02, 0.01, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};	
	const double tInit = 0.0;
	const double tEnd = 300.0;
	const int nt = 30001;
	
	// Define exogVar
	int nVarExogVar =2;
	int nSamplesVarExogVar[2] = {11, 3}; 
	double samplesExogVar[14] = {0,1,2,3,4,5,6,7,8,9,10,0, 50, 100};

	exogVarPeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	const int nV = 3;
	const int nIV = 5;

	GoodwinKeenRK4EV<double> keenModel(nV, nIV, tInit, tEnd, nt, &myExogVar); 
	// DEFINE CELLMAP
	int nCells1Dim = 151;
	int nCells[nV] = {nCells1Dim, nCells1Dim, nCells1Dim};
	double yLow[nV] = {-1.0, -1.0, -1.0};
	double yUp[nV] = {1.0, 1.0, 1.0};
	double byCells[nV];
	for (int it=0; it<nV; it++) {
		byCells[it] = (yUp[it] - yLow[it])/((double)(nCells1Dim-1))*(nCells[it]-1)/nCells[it];
	}

	double epsilon = 0.001;	
	double tol = 0.001;
	int nIterMax = 50;
	cellMapNR<double> myCellMapNR(&keenModel, parms, nCells, byCells, yLow, yUp, epsilon, tol, nIterMax);
	//NR<double> myNR(&keenModel, epsilon, tol, nIterMax);
	SCM<cellMapNR<double>> mySCMNR(myCellMapNR, false);
       // double y0[nV] = {0.836053, 0.968612, 0.0701911};
	//arma::Mat<double> J0(3, 3);	
	//myNR.Jacobian(y0, J0);
	//cout<<J0<<endl;
       // cout<<J0(1,1)<<endl;	


	// SIMULATION

	auto begin1=std::chrono::steady_clock::now();	
	mySCMNR.computeSCM();
	auto end1=std::chrono::steady_clock::now();
	
	auto timei=std::chrono::duration_cast<std::chrono::seconds>(end1-begin1).count();
	cout<<"Total computation time: "<<timei<<endl;
	// PLOT THE RESULTS
	std::vector<std::vector<int>> allAttractors = *mySCMNR.returnAttractors();
	std::cout<<"Raw SCM Algorithm: number of attractors: "<<allAttractors.size()<<endl;
	for (size_t it=0; it<allAttractors.size(); it++) {
		double tempY2[nV];
		myCellMapNR.cellToCoord(allAttractors[it][0], tempY2);
		cout<<"attractor number "<<it<<" cell: "<<allAttractors[it][0]<<" coord: "<<tempY2[0]<<" "<<tempY2[1]<<" "<<tempY2[2]<<endl;
	}
	
	// REFINE THE RESULTS
	mySCMNR.refineResults();

	// PLOT REFINED RESULTS
	allAttractors = *mySCMNR.returnAttractors();
	std::cout<<"SCM + Refine results: number of attractors: "<<allAttractors.size()<<endl;
	for (size_t it=0; it<allAttractors.size(); it++) {
		double tempY2[nV];
	      	myCellMapNR.cellToCoord(allAttractors[it][0], tempY2);
		cout<<"attractor number "<<it<<" cell: "<<allAttractors[it][0]<<" coord: "<<tempY2[0]<<" "<<tempY2[1]<<" "<<tempY2[2]<<endl;
	}

	return 0;
}	
