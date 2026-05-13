#ifndef CELLMAP_H
#define CELLMAP_H

#if UseRCPP==1
	#include "preprocRCPP_R.h"
#else
	#include "preproc.h"	
#endif

#include "ODE.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <stdio.h>
#include <string.h>

using namespace std;



template<typename T>
class cellMapBasin {
public:
	// CONSTRUCTOR
	cellMapBasin(ODE<T>* myModelIn, T* parmsIn,
		     int ntIn, T tEndIn,
		     const int* nCellsIn, const T* byCellsIn, 
		     const T* yLowIn, const T* yUpIn, const int tMultIn = 100) : 
			myModel(myModelIn), parms(parmsIn),
			nbCells(defineNbCells(nCellsIn)), 
			nCells(nCellsIn),
			byCells(byCellsIn),
			yLow(yLowIn),
			yUp(yUpIn),
			tMult(tMultIn),
		        sinkCell(defineNbCells(nCellsIn)-1) {
				myModel->changeNt(ntIn);
				myModel->changeTEnd(tEndIn);
				myModel->changeTInit(0.0);
			}

	// GET COORDINATES OF CELL cell
	void cellToCoord(int cell, T* y) {
		int multiplicator = 1;
		int coord[myModel->getNV()]={0};
		for (int it = 0; it<(myModel->getNV()-1); it++) {
			multiplicator*=nCells[it];
		}
		for (int it=myModel->getNV()-1; it>=0; it--) {
			while(cell >= multiplicator) {
				coord[it] += 1;
				cell-=multiplicator;
			}
			if (it>0) {
				multiplicator = multiplicator/nCells[it-1];
			}
		}
		for (int it=0; it<myModel->getNV(); it++) {
			y[it] = (coord[it]+0.5)*byCells[it] + yLow[it];
		}
	}

	// GET CELL NUMBER OF POINT y
	int coordToCell(T* y) {

		int cell = 0;
		int cellCoordinates[myModel->getNV()];
		int multiplicator = 1;

		for (int it=0; it<myModel->getNV(); it++) {
			if (y[it]>yUp[it] || y[it]<yLow[it] || isnan(y[it])) {
				return nbCells-1;
			}
			if (y[it]==yUp[it])
				y[it]-=byCells[it]/10.0f;
			if (y[it]==yLow[it])
				y[it]+=byCells[it]/10.0f;
			cellCoordinates[it] = (int) (floor((y[it] - yLow[it])/ byCells[it]));
		}

		for (int it=0; it<myModel->getNV(); it++) {
			cell += cellCoordinates[it]*multiplicator;
			multiplicator*=nCells[it];
		}
		return cell;
	}

	// GET IMAGE OF CELL cellIn
	int getCellImage(int cellIn) {
		T yInit[myModel->getNRowOut()], yEnd[myModel->getNRowOut()];
		int cellOut;
		cellToCoord(cellIn, yInit);
		myModel->solveLastPoint(yInit, parms, yEnd);
		cellOut = coordToCell(yEnd);
		//cout<<"init: "<<yInit[0]<<" "<<yInit[1]<<" "<<yInit[2]<<endl;
		//cout<<"end: "<<yEnd[0]<<" "<<yEnd[1]<<yEnd[2]<<endl;
		return cellOut;	
	}	

	// SAME AS ABOVE BUT WITH A LONGER TIME STEP
	int getCellImageRefined(int cellIn) {
		T y[myModel->getNRowOut()];
		int cellOut;
		cellToCoord(cellIn, y);
		for (int it=0; it<tMult; it++) {
			myModel->solveLastPoint(y, parms, y);
		}
		cellOut = coordToCell(y);

		return cellOut;	
	}
protected:
	ODE<T>* myModel;
	T* parms;
	const int nbCells;     // total number of cells
	const int* nCells;     // number of cells per myModel->getNV()ension
	const T* byCells;      // step length between two cells center
	const T* yLow;         // lower corner of cell
	const T* yUp;          // upper corner of cell
	const int tMult;      // How many times longer the simulation length should be for refined estimation of images
	const int sinkCell;

	int defineNbCells(const int* nCells) {
		int nbCells = 1;
		for (int it=0; it<myModel->getNV(); it++)
			nbCells*=nCells[it];
		nbCells+=1; // Add 1 for the sink cell
		return nbCells;
	}
};

#endif
