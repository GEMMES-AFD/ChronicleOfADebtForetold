#ifndef HCMSUBDIVISION_H
#define HCMSUBDIVISION_H

//#include "src/minDistRK4.h"

#if UseRCPP==1
	#include "preprocRCPP_R.h"
#else
	#include "preproc.h"	
#endif

#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <cstring>
#include <stdio.h>

using namespace std;

template<class cellMap, typename TInt>
class HCMSubdivision: public cellMap {
public:
	//CONSTRUCTOR
	HCMSubdivision(cellMap& aCellMap, int nSubdivisionMaxIn, int nSubdivisionGCMIn=0) : 
			cellMap{aCellMap}, 
			nSubdivisionMax(nSubdivisionMaxIn),
		        nSubdivisionGCM(nSubdivisionGCMIn) {
		
		// Initialize periodicCells and attractors (Initially one attractor: the sink cell)
		toExplore.resize(nbCells-1);
		nbCellsToExplore = nbCells;
		for (int it=0; it<nbCells-1; it++) {
			toExplore[it] = it;
		}
		cellMapImage = new TInt[nbCells];
		cellMapImage[sinkCell] = sinkCell;

		attractors.push_back(sinkCell);
		
		vector<TInt> temp(1);
		temp[0] = sinkCell;
		periodicCells.push_back(temp);
	}
	

	// DESTRUCTOR
	~HCMSubdivision() {
		delete[] cellMapImage;
	}

	
	// SEARCH FOR STATIC CELLS ONLY (I.E. CELLS THAT MAPS TO THEMSELVE)
	void computeHCMSubdivision() {
		int nSubDivision=0;
		while(nSubDivision<nSubdivisionMax && toExplore.size()>0) {
			if (nSubDivision<nSubdivisionGCM) {
				computeMapImageMultiple();
			} else {
				computeMapImage();
			}
			updateToExploreStatic();
			subdiviseCellMap();
			nSubDivision++;
		}
		computeMapImageRefined();
		for (int it=0; it<toExplore.size(); it++) {
			if(cellMapImage[it] == toExplore[it]) attractors.push_back(toExplore[it]);
		}
		// Remove sink cell
		attractors.erase(attractors.begin());
		if(toExplore.size()==0) {
			cout<<"No equilibrium found, at iteration "<<nSubDivision+1<<". Try to increase the initial number of cells, the density of the grid for GCM, or the number of iterations for which to use GCM."<<endl;
			cout<<"terminating."<<endl;
		}
	}

	

	//RETURN RESULTS
	// WARNING: DOES NOT PERFORM A DEEP COPY AND PROVIDES DIRECT ACCESS TO THE DATA STRUCTURE
	// DO NOT MODIFY IT
	// SUPPRESSING THE HCM CLASS OBJECT WILL DESTROY THE RESULTS	
	std::vector<TInt>* returnResults() {
		return &attractors;
	}	

protected:
	


	// INHERIT FUNCTIONS FROM PARENT CLASS
	using cellMap::cellMap;
	using cellMap::getCellImage; 
	using cellMap::getCellImageGCM; 
	using cellMap::getCellImageRefined;
	using cellMap::nbCells; //WARNING: THIS IS THE NUMBER OF CELLS IN THE FULL CELLMAP, NOT THE NUMBER OF CELLS THAT WILL BE EXPLORED!!
	using cellMap::sinkCell;
	using cellMap::subdiviseCellMap;
	using cellMap::getChildsId;	
	using cellMap::nChilds;

	// COMPUTE IMAGE OF ALL CELLS OF cellMap
	void computeMapImage() {

		#pragma omp parallel for if(UseParallel!=0) 
		for (int it=0; it<nbCellsToExplore-1; it++) {
			cellMapImage[it] = getCellImage(toExplore[it]);
		}
	}	
	
	// Same as Above, but with multiple starting points in each cell
	// Will work only for the research of periodic cells of order 1
	void computeMapImageMultiple() {

		#pragma omp parallel for if(UseParallel!=0) 
		for (int it=0; it<nbCellsToExplore-1; it++) {
			cellMapImage[it] = getCellImageGCM(toExplore[it]);
		}
	}
	// COMPUTE IMAGE OF ALL CELLS OF cellMap
	void computeMapImageRefined() {

		#pragma omp parallel for if(UseParallel!=0) 
		for (int it=0; it<nbCellsToExplore-1; it++) {
			cellMapImage[it] = getCellImageRefined(toExplore[it]);
		}
	}	

	// FIND PERIODIC CELLS OF ORDER 1 AND UPDATE LIST OF SELECTED CELLS
	 void updateToExploreStatic() {
		int comptSelectedCells = 0;
		TInt selectedCells[nbCellsToExplore-1];
		TInt childs[nChilds];
		for (int it=0; it<nbCellsToExplore-1; it++) {
			if(cellMapImage[it]==toExplore[it]) { // Select cells that map to themselves
				//double temp[3];
				//this->cellToCoord(toExplore[it], temp);
				//cout<<"parent: "<<toExplore[it]<<" "<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<endl;
				selectedCells[comptSelectedCells] = toExplore[it];
				comptSelectedCells++;
			}
		}
		//cout<<"number of selected cells: "<<comptSelectedCells<<" out of "<<nbCellsToExplore-1<<" cells"<<endl;
		//cout<<toExplore.size()<<endl;
		//selectedCells.erase(selectedCells.begin() + comptSelectedCells, selectedCells.end()); // Not needed
		toExplore.clear();
		toExplore.resize(comptSelectedCells*2);
		for (int it=0; it<comptSelectedCells; it++) {
			getChildsId(selectedCells[it], &toExplore[2*it]);
		}
		// Update cellMapImage
		delete[] cellMapImage;
		cellMapImage = new TInt[comptSelectedCells*2+1];
		cellMapImage[comptSelectedCells*2] = sinkCell;

		nbCellsToExplore = comptSelectedCells*2+1;
	 }
	

	TInt* cellMapImage;					// Image of all cells of the cellMap
	std::vector<TInt> toExplore;			// Id of cells to explore
	const int nSubdivisionMax;				// Max number of subdivisions of cellMap
	const int nSubdivisionGCM;				// Number of subdivision for which we use general and not simple cell mapping
	std::vector<TInt> attractors;			// ID of all attractors
	std::vector<std::vector<TInt>> periodicCells;	// Members of all attractors	
	TInt nbCellsToExplore;

};

#endif
