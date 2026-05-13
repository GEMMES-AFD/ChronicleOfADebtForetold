#ifndef CELLMAPNR_H
#define CELLMAPNR_H

#if UseRCPP==1
	#include "preprocRCPP_R.h"
#else
	#include "preproc.h"	
#endif

#include "newtonRaphson.h"
#include "ODE.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <stdio.h>
#include <string.h>

// WARNING: NOT PARALLELIZABLE BECAUSE THE SOLVE ALGORITHM FROM ARMA IS NOT THREAD SAFE 

using namespace std;

template<typename T>
class cellMapNR: public NR<T> {
public:
	// CONSTRUCTOR
	cellMapNR(ODE<T>* aModelIn, T* parmsIn,
		  const int* nCellsIn, const T* byCellsIn, 
		  const T* yLowIn, const T* yUpIn, 
		  const T epsilonIn, const T tolIn = 0.01, const int nIterMaxIn = 10) : 
			NR<T>{aModelIn, epsilonIn, tolIn, nIterMaxIn}, 
			parms(parmsIn),
			nbCells(defineNbCells(nCellsIn)), 
			nCells(nCellsIn),
			byCells(byCellsIn),
			yLow(yLowIn),
			yUp(yUpIn),
		        sinkCell(defineNbCells(nCellsIn)-1) {}

	// INHERIT FUNCTIONS FROM MODEL<T>
	using NR<T>::NR;
	using NR<T>::NR1Iter; 
	using NR<T>::NRMethod; 


	// GET COORDINATES OF CELL cell
	void cellToCoord(int cell,
		         T* y) {
		int multiplicator = 1;
		int coord[this->getNV()]={0};
		for (int it = 0; it<(this->getNV()-1); it++) {
			multiplicator*=nCells[it];
		}
		for (int it=this->getNV()-1; it>=0; it--) {
			while(cell >= multiplicator) {
				coord[it] += 1;
				cell-=multiplicator;
			}
			if (it>0) {
				multiplicator = multiplicator/nCells[it-1];
			}
		}
		for (int it=0; it<this->getNV(); it++) {
			y[it] = (coord[it]+0.5)*byCells[it] + yLow[it];
		}
	}

	// GET CELL NUMBER OF POINT y
	int coordToCell(T* y) {

		int cell = 0;
		int cellCoordinates[this->getNV()];
		int multiplicator = 1;

		for (int it=0; it<this->getNV(); it++) {
			if (y[it]>yUp[it] || y[it]<yLow[it] || isnan(y[it])) {
				return nbCells-1;
			}
			if (y[it]==yUp[it])
				y[it]-=byCells[it]/10.0f;
			if (y[it]==yLow[it])
				y[it]+=byCells[it]/10.0f;
			cellCoordinates[it] = (int) (floor((y[it] - yLow[it])/ byCells[it]));
		}

		for (int it=0; it<this->getNV(); it++) {
			cell += cellCoordinates[it]*multiplicator;
			multiplicator*=nCells[it];
		}
		return cell;
	}

	// GET IMAGE OF CELL cellIn
	int getCellImage(int cellIn) {
		T y[this->getNRowOut()];
		int cellOut;
		cellToCoord(cellIn, y);
		NR1Iter(this->getTInit(), y, parms);
		cellOut = coordToCell(y);
		//cout<<"init: "<<yInit[0]<<" "<<yInit[1]<<" "<<yInit[2]<<endl;
		//cout<<"end: "<<yEnd[0]<<" "<<yEnd[1]<<yEnd[2]<<endl;
		return cellOut;	
	}	
	
	// Perform full Newton Raphson algorithm to find exact coordinates of zeroes of the system.
	int getCellImageRefined(int cellIn) {
		T y[this->getNRowOut()];
		int cellOut;
		cellToCoord(cellIn, y);
		NRMethod(this->getTInit(), y, parms);
		cellOut = coordToCell(y);

		return cellOut;
	}

protected:
	T* parms;
	const int nbCells;     // total number of cells
	const int* nCells;     // number of cells per dimension
	const T* byCells;      // step length between two cells center
	const T* yLow;         // lower corner of cell
	const T* yUp;          // upper corner of cell
	const int sinkCell;

	int defineNbCells(const int* nCells) {
		int nbCells = 1;
		for (int it=0; it<this->getNV(); it++)
			nbCells*=nCells[it];
		nbCells+=1; // Add 1 for the sink cell
		return nbCells;
	}
};

#endif
