#ifndef modelR_h
#define modelR_h

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <stdio.h>
#include <string.h>

#if UseRCPP==1
        #include "preprocRCPP_R.h"
#else
        #include "preproc.h"    
#endif

#include "ODE.h"
#include "exogenousVariables.h"

using namespace std;

// GOODWIN-KEEN BUT WITH TWO EXOGENOUS VARIABLES (THIS IS JUST TO ILLUSTRATE THE USE OF INTERMEDIATE VARIABLES)
template<typename T>
class modelR : virtual public ODE<T> {
public:
	@ADDEventTime

	@ADDEventVar

	// GOODWIN-KEEN MODEL DEFINITION
	@ADDFunc	
};

#endif
