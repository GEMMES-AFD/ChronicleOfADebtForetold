#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <cblas.h>
#define ARMA_WARN_LEVEL 1
#include <armadillo>

#include "src/exogenousVariables.h"

double truePoly(double t, int i) {
	if (i==0 && t<=1) return -0.3057*pow(t, 3.0) + 3.3036*t + 21;
	if (i==0 && t<=2) return -1.4821*pow(t, 3.0) + 3.5357*pow(t, 2) - 0.23214*t + 22.179;
	if (i==0 && t<=3) return 3.2321*pow(t, 3.0) - 24.750*pow(t, 2) + 56.339*t - 15.536;
	if (i==0 && t<=4) return -1.4464*pow(t, 3.0) + 17.357*pow(t, 2) -69.982*t + 110.79;
	if (i==1) return t*2.0;
	return -1.0;
}

void compareResult(double t, int i, exogVarCubicSplinePeriodic<double>& myExogVar) {
	cout<<"i = "<<i<<" t = "<<t<<" true sol: "<<truePoly(t, i)<<" estimated poly: "<<myExogVar.getValue(t, i)<<endl;
}	

int main() {
	const double tInit = 0.0;
	const double tEnd = 4;

	// Define exogVar
	double samplesExogVar[14] = {21,24,24,18,16, 0, 1, 2, 3, 4, 5, 6, 7, 8};
	int nVarExogVar = 2;
	int nSamplesVarExogVar[nVarExogVar] = {5, 9}; 

	exogVarCubicSplinePeriodic<double> myExogVar(tInit, tEnd, samplesExogVar, nVarExogVar, nSamplesVarExogVar);
	int compt=0;
	for(int j=0; j<myExogVar.nVar; j++) {
		for (int i=0; i<myExogVar.nSamplesVar[j]-1; i++) {
			cout<<"a: "<<myExogVar.a[compt+i]<<" b "<<myExogVar.b[compt+i]<<" c "<<myExogVar.c[compt+i]<<" d "<<myExogVar.d[compt+i]<<endl;
		}
		compt+=myExogVar.nSamplesVar[j]-1;
	}
	compareResult(tInit, 0, myExogVar);
	compareResult(tEnd, 0, myExogVar);
	compareResult(0.99, 0, myExogVar);
	compareResult(1.9, 0, myExogVar);

	compareResult(tInit, 1, myExogVar);
	compareResult(tEnd, 1, myExogVar);
	compareResult(0.99, 1, myExogVar);
	compareResult(1.9, 1, myExogVar);

	return 0;
}
