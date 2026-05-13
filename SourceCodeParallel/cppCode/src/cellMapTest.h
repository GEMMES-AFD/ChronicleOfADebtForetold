#ifndef CELLMAPTEST_H
#define CELLMAPTEST_H

#if UseRCPP==1
	#include "preprocRCPP_R.h"
#else
	#include "preproc.h"	
#endif


#include "RK4.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <stdio.h>
#include <string.h>

using namespace std;



template<typename T>
class cellMapTest {
public:
	// CONSTRUCTOR
	cellMapTest(const int* nCellsIn) : 
			nbCells(defineNbCells(nCellsIn)), 
			nCells(nCellsIn),
		        sinkCell(defineNbCells(nCellsIn)-1) {
			}

	// GET IMAGE OF CELL cellIn
	int getCellImage(int cellIn) {
		if(cellIn<10 && cellIn>1) return cellIn-1;
		if(cellIn<=1)return 9;

		if (cellIn<50 && cellIn>40) return cellIn+1;
		if (cellIn==50) return 5;

		return cellIn+1;	
	}	
	// SAME AS ABOVE BUT WITH A LONGER TIME STEP
	int getCellImageRefined(int cellIn) {
		return getCellImage(cellIn);	
	}
protected:
	const int nbCells;     // total number of cells
	const int* nCells;     // number of cells per dimension
	const int sinkCell;
	int defineNbCells(const int* nCells) {
		int nbCells = 1;
		for (int it=0; it<dim; it++)
			nbCells*=nCells[it];
		nbCells+=1; // Add 1 for the sink cell
		return nbCells;
	}
};

#endif
