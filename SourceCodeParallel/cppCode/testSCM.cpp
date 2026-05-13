#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

#include "src/preproc.h"
#include "src/SCM.h"
#include "src/cellMapBasin.h"
#include "src/goodwinKeenRK4.h"

using namespace std;


int main() {
	// INIT MAIN PARAMETERS
	const int parmsLength = 10;
	double parms[parmsLength] = {0.025, 0.02, 0.01, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};	
	const double tInit = 0.0;
	const double tEnd = 10.0;
	const int nt = 301;
	
	// Define exogVar
	int nVarExogVar = 2;
	int nSamplesVarExogVar[2] = {11, 3}; 
	double samplesExogVar[14] = {0,1,2,3,4,5,6,7,8,9,10,0, 50, 100};

	exogVarPeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	const int nV = 3;
	const int nIV = 5;

	GoodwinKeenRK4EV<double> keenModel(nV, nIV, tInit, tEnd, nt, &myExogVar); 
	
	// DEFINE CELLMAP
	int nCells1Dim = 51;
	int nCells[nV] = {nCells1Dim, nCells1Dim, nCells1Dim};
	double yLow[nV] = {0.0, 0.0, 0.0};
	double yUp[nV] = {2.0, 1.0, 10.0};
	double byCells[nV];
	for (int it=0; it<nV; it++) {
		byCells[it] = (yUp[it] - yLow[it])/((double)(nCells1Dim-1))*(nCells[it]-1)/nCells[it];
	}
	int ntCellMap = 500;
	int ntMult = 50000/ntCellMap;
	cellMapBasin<double> myCellMapBasin(&keenModel, parms, ntCellMap, tEnd, nCells, byCells, yLow, yUp, ntMult);

	SCM<cellMapBasin<double>> mySCM(myCellMapBasin);

	// SIMULATION
	mySCM.computeSCM();

	// PLOT THE RESULTS
	std::vector<std::vector<int>> allAttractors = *mySCM.returnAttractors();
	std::vector<std::vector<int>> allDomainsOfAttraction = *mySCM.returnDomainsOfAttraction();
	std::cout<<"Raw SCM Algorithm: number of attractors: "<<allAttractors.size()<<endl;
	for (size_t it=0; it<allAttractors.size(); it++) {
		double tempY2[nV];
	      	myCellMapBasin.cellToCoord(allAttractors[it][0], tempY2);
		cout<<"attractor number "<<it<<"cell: "<<allAttractors[it][0]<<"coord: "<<tempY2[0]<<" "<<tempY2[1]<<" "<<tempY2[2]<<endl;
		cout<<"number of cells in domain of attraction: "<<allDomainsOfAttraction[it].size()<<endl;
	}

	// REFINE THE RESULTS
	mySCM.refineResults();

	// PLOT REFINED RESULTS
	allAttractors = *mySCM.returnAttractors();
	allDomainsOfAttraction = *mySCM.returnDomainsOfAttraction();
	std::cout<<"SCM + Refine results: number of attractors: "<<allAttractors.size()<<endl;
	for (size_t it=0; it<allAttractors.size(); it++) {
		double tempY2[nV];
	      	myCellMapBasin.cellToCoord(allAttractors[it][0], tempY2);
		cout<<"attractor number "<<it<<" cell: "<<allAttractors[it][0]<<" coord: "<<tempY2[0]<<" "<<tempY2[1]<<" "<<tempY2[2]<<endl;
		cout<<"number of cells in domain of attraction: "<<allDomainsOfAttraction[it].size()<<endl;
	}


	return 0;
}	
