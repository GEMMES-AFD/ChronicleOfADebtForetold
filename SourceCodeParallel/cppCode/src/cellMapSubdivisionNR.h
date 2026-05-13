#ifndef CELLMAPSUBDIVISIONNR_H
#define CELLMAPSUBDIVISIONNR_H

#if UseRCPP==1
	#include "preprocRCPP_R.h"
#else
	#include "preproc.h"	
#endif

#include "newtonRaphson.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <stdio.h>
#include <string.h>

using namespace std;


template<typename T, typename TInt>
class cellMapSubdivisionNR: public NR<T> {
public:
	// CONSTRUCTOR
	cellMapSubdivisionNR(ODE<T>* aModelIn, T* parmsIn,
			     int* nCellsIn, T* byCellsIn,
			     const T* yLowIn, const T* yUpIn, 
			     const T epsilonIn, const T tolIn = 0.01, 
			     const int nIterMaxIn = 10,
			     T* gridCellIn = nullptr, 
			     int nEltsGridCellIn = 0) : 
			NR<T>{aModelIn, epsilonIn, tolIn, nIterMaxIn}, 
			parms(parmsIn),
			nbCells(defineNbCells(nCellsIn)), 
			byCells(byCellsIn),
			yLow(yLowIn),
			yUp(yUpIn),
		        sinkCell(defineNbCells(nCellsIn)-1),
		        dimToSubdivise(0), 
		       	nEltsGridCell(nEltsGridCellIn) {
		gridCell = new T[nEltsGridCellIn*this->getNV()];
		for (TInt it=0; it<nEltsGridCellIn*this->getNV(); it++) {
			gridCell[it] = gridCellIn[it];	
		}
		nCells = new int[this->getNV()];
		for (int it=0; it<this->getNV(); it++) nCells[it]=nCellsIn[it];
	}

	// COPY CONSTRUCTOR
	cellMapSubdivisionNR(cellMapSubdivisionNR<T, TInt> const& other) : 
			NR<T>{other},
			parms(other.parms),
			nbCells(other.nbCells), 
			byCells(other.byCells),
			yLow(other.yLow),
			yUp(other.yUp),
		        sinkCell(other.sinkCell),
		        dimToSubdivise(other.dimToSubdivise), 
		       	nEltsGridCell(other.nEltsGridCell) {
		gridCell = new T[other.nEltsGridCell*this->getNV()];
		for (TInt it=0; it<other.nEltsGridCell*this->getNV(); it++) {
			gridCell[it] = other.gridCell[it];	
		}
		nCells = new int[this->getNV()];
		for (int it=0; it<this->getNV(); it++) nCells[it]=other.nCells[it];

	} 
	
	// DESTRUCTOR
	~cellMapSubdivisionNR() {
		delete[] gridCell;
		delete[] nCells;
	}



	// INHERIT FUNCTIONS
	using NR<T>::NR;
	using NR<T>::NR1Iter; 
	using NR<T>::NRMethod; 


	//************************************************************************//
	//** There are three ways to define a cell:                             **//
	//** (1) Its scalar id: int cell                                        **//
	//** (2) Its integer coordinates in the cell space int[dim] cellCoord   **//
	//** (3) The coordinates of the cell center in the phase space T[dim] y **//
	//** The functions below are use to perform conversion from             **//
	//** one type to another:                                               **//
	//** void cellScalarToCellVec: (1) -> (2)                               **//
	//** int cellVecToCellScalar:  (2) -> (1)                               **//
	//** void cellToCoord:         (1) -> (3)                               **//
	//** int coordToCell:          (3) -> (1)                               **//
	//************************************************************************//

	void cellScalarToCellVec(TInt cell, int* cellCoord) {
		TInt multiplicator = 1;
		for (int it = 0; it<(this->getNV()-1); it++) {
			multiplicator*=nCells[it];
		}
		for (int it=this->getNV()-1; it>=0; it--) {
			while(cell >= multiplicator) {
				cellCoord[it] += 1;
				cell-=multiplicator;
			}
			if (it>0) {
				multiplicator = multiplicator/nCells[it-1];
			}
		}
	}

	TInt cellVecToCellScalar(int* cellCoord) {
		TInt cell = 0;
		TInt multiplicator = 1;
		for (int it=0; it<this->getNV(); it++) {
			cell += cellCoord[it]*multiplicator;
			multiplicator*=nCells[it];
		}
		return cell;
	}

	// GET COORDINATES OF CELL cell
	void cellToCoord(TInt cell,
		      T* y) {
		int cellCoord[this->getNV()]={0};
		cellScalarToCellVec(cell, cellCoord);
		for (int it=0; it<this->getNV(); it++) {
			y[it] = (cellCoord[it]+0.5)*byCells[it] + yLow[it];
		}
	}

	// GET CELL NUMBER OF POINT y
	TInt coordToCell(T* y) {
		int cellCoordinates[this->getNV()];

		for (int it=0; it<this->getNV(); it++) {
			if (y[it]>yUp[it] || y[it]<yLow[it] || isnan(y[it])) {
				return sinkCell;
			}
			if (y[it]==yUp[it]) {
				y[it]-=byCells[it]/10.0f;
			}
			if (y[it]==yLow[it])
				y[it]+=byCells[it]/10.0f;
			cellCoordinates[it] = (int) (floor((y[it] - yLow[it])/ byCells[it]));
		}
			return cellVecToCellScalar(cellCoordinates);
	}
	
	// Get cell id after performing one more subdivision
	TInt cellVecToCellScalarNextGrid(int* cellCoord) {
		TInt cell = 0;
		TInt multiplicator = 1;
		for (int it=0; it<this->getNV(); it++) {
			cell += cellCoord[it]*multiplicator;
			multiplicator*=nCells[it];
			if (it==dimToSubdivise) multiplicator*=2;
		}
		return cell;
	}

	// GET IMAGE OF CELL cellIn
	TInt getCellImage(TInt cellIn) {
		T y[this->getNRowOut()];
		TInt cellOut;
		cellToCoord(cellIn, y);
		NR1Iter(this->getTInit(), y, parms);
		cellOut = coordToCell(y);
		return cellOut;	
	}

	// Same as above, but with multiple points in each cell
	// WARNING // Only finds periodic cells of period 1
	// Returns current cell if at least one point remains in current cell
	// Otherwise, returns the sink cell
	TInt getCellImageGCM(TInt cellIn) {
		T y0[this->getNRowOut()];
		cellToCoord(cellIn, y0);
		for (int it=0; it<this->getNV(); it++) y0[it]-=byCells[it]/2;
		T y[this->getNRowOut()];
		TInt cellOut;
		for (TInt it=0; it<nEltsGridCell; it++) {
			for (int it2=0; it2<this->getNV(); it2++) {
				y[it2] = y0[it2] + (gridCell[it*this->getNV() + it2])*byCells[it2];
			}
			NR1Iter(this->getTInit(), y, parms);
			cellOut = coordToCell(y);
			if (cellOut==cellIn) return cellIn;
		}	
		return sinkCell;
	}


	// Same as getCellImage, but performs multiples iterations of the Newton-Raphson method
	TInt getCellImageRefined(TInt cellIn) {
		T y[this->getNRowOut()];
		TInt cellOut;
		cellToCoord(cellIn, y);
		NRMethod(this->getTInit(), y, parms);
		cellOut = coordToCell(y);

		return cellOut;
	}
protected:
	
	// Double the number of cells by cutting each cell in two along the dimension dimToSubdivise
	void subdiviseCellMap() {
		nCells[dimToSubdivise]*=2;
		byCells[dimToSubdivise]/=2;
		nbCells = defineNbCells(nCells);
		sinkCell = nbCells-1;

		// Update dimension along which to perform (next) subdivision;
		dimToSubdivise++;
		if(dimToSubdivise>=this->getNV()) dimToSubdivise=0;
	}	

	// RETURNS ID OF cell'S CHILDRENS AFTER SUBDIVISION	
	void getChildsId(TInt cell, TInt* childs) {
		int parentCoord[this->getNV()]={0};
		cellScalarToCellVec(cell, parentCoord);
		parentCoord[dimToSubdivise] = parentCoord[dimToSubdivise]*2;
		for (int it=0; it<nChilds; it++) {
			childs[it] = cellVecToCellScalarNextGrid(parentCoord);
		        parentCoord[dimToSubdivise]++;	
		}	
	}

	// Return the total number of cells in the cellMap
	TInt defineNbCells(const int* nCells) {
		TInt out = 1;
		for (int it=0; it<this->getNV(); it++)
			out*=nCells[it];
		out+=1; // Add 1 for the sink cell
		return out;
	}
	T* parms;
	TInt nbCells;   // total number of cells
	int* nCells;             // number of cells per dimension
	T* byCells;              // step length between two cells center
	const T* yLow;           // lower corner of cell
	const T* yUp;            // upper corner of cell
	TInt sinkCell;	         // Coordinates of the sinkCell
	int dimToSubdivise;	 // Dimension along which to perform subdivision
	int nChilds = 2;	 // Number of child per cell during subdivision
	int nEltsGridCell;	 // Number of points per cell for general mapping
	T* gridCell;		 // Coordinates of points in each cell (i.e. in a [0; 1]^this->getNV() square)
};

#endif
