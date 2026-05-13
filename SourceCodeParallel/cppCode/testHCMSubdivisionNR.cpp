#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define ARMA_WARN_LEVEL 1
#include <armadillo>


#include "src/preproc.h"
#include "src/goodwinKeenRK4.h"
#include "src/HCMSubdivision.h"
#include "src/cellMapSubdivisionNR.h"

using namespace std;
typedef long long unsigned int TInt;
typedef double T;
void cellScalarToCellVec(int cell, int nCells, int* cellCoord, int nV) {
	int multiplicator = pow(nCells, nV-1);
	for (int it=nV-1; it>=0; it--) {
		cellCoord[it] = 0;
		while(cell >= multiplicator) {
			cellCoord[it] += 1;
			cell-=multiplicator;
		}
		if (it>0) {
			multiplicator = multiplicator/nCells;
		}
	}
}

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

	// DEFINE CELLMAP
	int nCells1Dim = 5;
	int nCells[nV];
	T yLow[nV];
	T yUp[nV];
	T byCells[nV];
	for (int it=0; it<nV; it++) {
		nCells[it] = nCells1Dim;
		yLow[it] = 0;
		yUp[it] = 1.0;
		byCells[it] = (yUp[it] - yLow[it])/((T)(nCells1Dim-1))*(nCells[it]-1)/nCells[it];
	}

	// Parms for Newton Raphson
	T epsilon = 0.00001;	        // step to estimate Jacobian matrix numerically
	T tol = 0.0001;	                // Tolerance for convergence criterion
	int nIterMax = 2000;            // Max number of iterations

	// Parms for subdivision and hybrid cell mapping
	int nSubdivisionMax = nV*4;               // Maximum number of subdivisions
	int nSubdivisionGCM = nSubdivisionMax-nV;// Number of subdivisions for which generall cell mapping is used (after the algorithm uses simple cell mapping)
	int nEltsPerDimGridCell = 5;		
	int nEltsGridCell = pow(nEltsPerDimGridCell, nV); // Number of points per cell for GCM
	T gridCell[nEltsGridCell*nV];                // Grid for each cell (defined on a [0, 1]^nV cube)
	int myIterators[nV]={0};
	for (int it=0; it<nEltsGridCell; it++) {
		cellScalarToCellVec(it, nEltsPerDimGridCell, myIterators, nV);
		for (int it2=0; it2<nV; it2++) {
			gridCell[it*nV + it2] = ((T) (myIterators[it2]+1))/((T) (nEltsPerDimGridCell));
		}
	}
	cellMapSubdivisionNR<T, TInt> myCellMapSubNR(&keenModel, parms, 
						     nCells, byCells, 
						     yLow, yUp, epsilon, tol, nIterMax, 
						     gridCell, nEltsGridCell);

	HCMSubdivision<cellMapSubdivisionNR<T, TInt>, TInt> myHCMSubNR(myCellMapSubNR, nSubdivisionMax, nSubdivisionGCM);
		
	// SIMULATION
	myHCMSubNR.computeHCMSubdivision();

	// PLOT THE RESULTS
	std::vector<TInt> allAttractors = *myHCMSubNR.returnResults();
	std::cout<<"Raw SCM Algorithm: number of equilibria: "<<allAttractors.size()<<endl;
	for (size_t it=0; it<allAttractors.size(); it++) {
		T tempY2[nV];
	        myHCMSubNR.cellToCoord(allAttractors[it], tempY2);
		myCellMapSubNR.NRMethod(tInit, tempY2, parms); // get "exact" equilibrium using Newton-Raphson
		cout<<"equilibrium number "<<it<<" cell: "<<allAttractors[it]<<" coord: ";
		for (int it2=0; it2<nV; it2++) cout<<tempY2[it2]<<" ";
		cout<<endl;
	}
	
	//// REFINE THE RESULTS
	//mySCMNR.refineResults();

	//// PLOT REFINED RESULTS
	//allAttractors = *mySCMNR.returnAttractors();
	//std::cout<<"SCM + Refine results: number of attractors: "<<allAttractors.size()<<endl;
	//for (size_t it=0; it<allAttractors.size(); it++) {
		//T tempY2[nV];
		      //myCellMapNR.cellToCoord(allAttractors[it][0], tempY2);
		//cout<<"attractor number "<<it<<" cell: "<<allAttractors[it][0]<<" coord: "<<tempY2[0]<<" "<<tempY2[1]<<" "<<tempY2[2]<<endl;
	//}

	return 0;
}	
