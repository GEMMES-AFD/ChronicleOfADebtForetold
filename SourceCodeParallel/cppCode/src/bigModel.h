#ifndef longModel_h
#define longModel_h

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

using namespace std;

template<typename T>
class bigModel : virtual public ODE<T> {
public:
	void makeEventTime(const T t, T* parms, T* y, T* x) override {} 

	void makeEventVar(const T t, T* parms, T* y, T* x) override {}

	void Func(const T t, const T* y, const T* parms, T* ydot, T* x) override {
		for (int it=0; it<100; it++) ydot[it] = sin(t);
		}
};

#endif
