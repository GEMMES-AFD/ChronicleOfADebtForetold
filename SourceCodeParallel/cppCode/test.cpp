#define UseRCPP 1 		// true if the code is called from R(via RCPP), false if it is directly called from c++


#include "src/preprocRCPP_R.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <Rcpp.h>

#define ARMA_WARN_LEVEL 0
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::plugins("cpp14")]]

using namespace std;

// [[Rcpp::export]]
double useless() {
	return 1.0;
	}
