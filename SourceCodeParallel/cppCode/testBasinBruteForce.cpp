#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <armadillo>

#include "src/preproc.h"
#include "src/basinBruteForce.h"
#include "src/goodwinKeenRK4.h"
using namespace std;


int main() {
	// INIT MAIN PARAMETERS
	const int parmsLength = 10;
	double parms[parmsLength] = {0.025, 0.02, 0.01, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};	
	const double tInit = 0;
	const double tEnd = 300;
	const int nt = 30001;
	
	const int nV = 3;
	const int nIV = 5;
	
	// Define exogVar
	int nVarExogVar = 2;
	int nSamplesVarExogVar[2] = {11, 3}; 
	double samplesExogVar[14] = {0,1,2,3,4,5,6,7,8,9,10,0, 50, 100};

	exogVarPeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	
	 


	// INIT PARAMETERS FOR BASIN OF ATTRACTION
	double yeq[3] = {0.83605402, 0.96861321, 0.07019744};
	double tol= 0.1;
	double boundsSD[2*nV] = {0.4,1.0, 0.4, 1.0, 0.0, 1.0};
	bool useDist=true;
	bool useForConv[nV] = {true, true, true};
	// DEFINE GRID TO EXPLORE
	int gridSize1Dim = 30;
	int gridSize = pow(gridSize1Dim, nV);
	cout<<"gridSize: "<<gridSize<<endl;
	double boundsLower[nV] = {0.0, 0.0, 0.0};
	double boundsUpper[nV] = {2.0, 1.0, 10.0};
	double increment[gridSize1Dim];
	for (int it=0; it<nV; it++) {
		increment[it] = (boundsUpper[it] - boundsLower[it])/((double)(gridSize1Dim-1));
	}

	double grid[gridSize*nV];
	int it1, it2, it3;
	int** myIterators = (int**)malloc(sizeof(int*)*nV);
	myIterators[0]=&it1;
	myIterators[1]=&it2;
	myIterators[2]=&it3;
	// ugly but it works !
	for (it1=0; it1<gridSize1Dim; it1++) {
		for (it2=0; it2<gridSize1Dim; it2++) {
			for (it3=0; it3<gridSize1Dim; it3++) {
				for (int it=0; it<nV; it++) {
					grid[it1*gridSize1Dim*gridSize1Dim*nV+it2*gridSize1Dim*nV+it3*nV+it] = increment[it]*(*myIterators[it]);
				}
			}
		}
	}
	// INIT DATA STRUCTS FOR BASIN COMPUTATION
	GoodwinKeenRK4EV<double> keenModel(nV, nIV, tInit, tEnd, nt, &myExogVar); 
	
	convModel<double> myConv= {useDist, useForConv, yeq, tol, boundsSD};
	gridModel<double> myGrid= {grid, gridSize};
	int out[gridSize];
	int comptPointsInBasin = 0;

	basin<double> myBasin(&keenModel, parms, myGrid, myConv);
	myBasin.computeBasin(out);
	
	for (int it=0; it<gridSize; it++) {
		if (out[it]==1) comptPointsInBasin++;
	}
	cout<<"number of points in bassin: "<<comptPointsInBasin<<endl;
	// FREE POINTERS
	free (myIterators);

	return 0;
}	
