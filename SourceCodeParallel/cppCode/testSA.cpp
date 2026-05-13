#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <armadillo>

#include "src/goodwinKeenRK4.h"
#include "src/preproc.h"
#include "src/sensitivityAnalysis.h"

using namespace std;

int main() {
	// INIT MAIN PARAMETERS
	const int parmsLength = 10;
	double parms[parmsLength] = {0.025, 0.02, 0.01, 3.0, 0.03, -0.0065, -5.0, 20.0, 0.04/(1.0-pow(0.04, 2.0)), pow(0.04, 3.0)/(1.0-pow(0.04, 2.0))};	
	const double tInit = 0.0;
	const double tEnd = 300.0;
	const int nt = 3001;
	
	// Define exogVar
	int nVarExogVar = 2;
	int nSamplesVarExogVar[2] = {11, 3}; 
	double samplesExogVar[14] = {0,1,2,3,4,5,6,7,8,9,10,0, 50, 100};

	exogVarPeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	const int nV = 3;
	const int nIV = 5;

	/*********************************/
	/*********************************/
	/* CODE FOR SENSITIVITY ANALYSIS */
	/*********************************/
	/*********************************/
	
	const int nParmsSetSA = 10;
	const int nParmsSA = 3;
	int parmsPosSA[nParmsSA] = {0, 3, 5};
	double* allParmsSA= (double*) malloc(sizeof(double)*nParmsSA*nParmsSetSA);
	for (int it=0; it<nParmsSetSA; it++) {
		for (int it2=0; it2<nParmsSA; it2++) {
			allParmsSA[it*nParmsSA + it2] = parms[parmsPosSA[it2]];
		}
		allParmsSA[it*nParmsSA + 0] = (double) 0.015 + (0.035-0.015)*it/nParmsSetSA;
	}
	double y0[dim] = {0.8, 0.9, 0.1};


	GoodwinKeenRK4EV<double> keenModel(nV, nIV, tInit, tEnd, nt, &myExogVar); 
	SA<double> mySA(&keenModel, parms, parmsLength, nParmsSetSA, nParmsSA, parmsPosSA, allParmsSA); 

	int nRowOut = keenModel.getNRowOut();	
	
	// RECOVER LAST POINT ONLY //	
	double* outSA = (double*) malloc(sizeof(double)*nRowOut*nParmsSetSA);

	mySA.sensitivityAnalysisLastPoint(y0, outSA);

	for (int it=0; it<nParmsSetSA; it++) {
		cout <<"Parameters: ";
			for (int it2=0; it2<nParmsSA; it2++) {       
				cout << allParmsSA[it*nParmsSA + it2] <<" ";
			}
			cout << "end position: ";
			for (int it2=0; it2<nRowOut; it2++) {	
				cout << outSA[it*nRowOut+ it2]<<" ";
			}
		cout<<endl;	
	}
	
	free(outSA);

	// RECOVER FULL TRAJECTORY //
	
	double* outSAFull= (double*) malloc(sizeof(double)*nRowOut*nt*nParmsSetSA);
	mySA.sensitivityAnalysisFullTrajectory(y0, outSAFull);
	for (int it=0; it<nParmsSetSA; it++) {
		cout <<"Parameters: ";
			for (int it2=0; it2<nParmsSA; it2++) {      
				cout << allParmsSA[it*nParmsSA + it2] <<" ";

			}
			cout<<endl;
			for (int it2=0; it2<nRowOut; it2++) {
				cout << "Trajectory for variable "<<it2<<": ";
				for (int it3=0; it3<(int)nt/1000+1; it3++) {
					cout << outSAFull[it*nRowOut*nt + it3*1000*nRowOut+ it2]<<" ";
				}
				cout<<endl;
			}
		cout<<endl;	
	}
	free(outSAFull);
	free(allParmsSA);

	return 0;
}
