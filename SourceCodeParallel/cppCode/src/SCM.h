#ifndef SCM_H
#define SCM_H

#if UseRCPP==1
	#include "preprocRCPP_R.h"
#else
	#include "preproc.h"	
#endif
#include <chrono>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <cstring>
#include <stdio.h>

#include "cellMapBasin.h"
#include "cellMapNR.h"

using namespace std;

template<class cellMap>
class SCM: public cellMap {
public:
	//CONSTRUCTOR
	SCM(cellMap& aCellMap, bool computeDomainOfAttractionIn=true) : 
			cellMap{aCellMap}, 
			computeDomainOfAttraction(computeDomainOfAttractionIn) {
		// Initialize periodicCells and attractors (only one attractor: the sink cell)
		cellMapImage = new int[nbCells];
		attractors.push_back(sinkCell);
		std::vector<int> temp(1);
		temp[0] = sinkCell;
		periodicCells.push_back(temp);
	}
	

	// DESTRUCTOR
	~SCM() {
		delete[] cellMapImage;
	}

	
	//MAIN FUNCTION: FINDS ATTRACTORS AND THEIR DOMAIN OF ATTRACTION
	void computeSCM() {
		computeMapImage();
		findAttractors();
		if(computeDomainOfAttraction) findDomainsOfAttraction();		
		// Remove sink cell
		attractors.erase(attractors.begin());
		periodicCells.erase(periodicCells.begin());
		if(computeDomainOfAttraction==true) domainsOfAttraction.erase(domainsOfAttraction.begin());
	}	

	//SAME AS ABOVE, BUT REFINE RESULTS TO REMOVE/MERGE "SPURIOUS" ATTRACTORS
	void computeSCM(bool refine) {
		computeSCM();
		if (refine==true) { 
			refineResults();
		}
	}
	

	// FIND AND MERGE/ELIMINATE "SPURIOUS" ATTRACTORS
	void refineResults() {
		std::vector<bool> mergedAttractors(attractors.size());
		std::vector<int> refinedImg(attractors.size());

		#pragma omp parallel for default(shared) if(UseParallel!=0)
		for (int it=0; it<attractors.size(); it++) {
						mergedAttractors[it] = false;
			refinedImg[it] = cellMapImage[getCellImageRefined(attractors[it])];
			//int refinedImg[it] = getCellImageRefined(attractors[it]);
		}	
		

		for(int  it=0; it<attractors.size(); it++) {
			//if(refinedImg[it]!=sinkCell) cout<<"it: "<<it<<" attractorCell: "<<attractors[it]<<" imageRefined: "<<refinedImg[it]<<endl;
			for (int it2=0; it2<attractors.size(); it2++) {
				 if(isInAttractor(refinedImg[it], it2) && it!=it2) {
				 		mergeAttractors(it2, it);
					mergedAttractors[it]=true;
				}
			}	
			if(refinedImg[it]==sinkCell) mergedAttractors[it] = true; // set attractors pointing to sinkCell as merged so they are erased below
		}
		// Delete merged attractors
		for (int it=mergedAttractors.size()-1; it>=0; it--) {
			if (mergedAttractors[it]==true) {
				periodicCells.erase(periodicCells.begin()+it);
				attractors.erase(attractors.begin()+it);
				if (computeDomainOfAttraction==true) domainsOfAttraction.erase(domainsOfAttraction.begin()+it);
			}
		}
	}
	
	
	//RETURN RESULTS
	// WARNING: DOES NOT PERFORM A DEEP COPY AND PROVIDES DIRECT ACCESS TO THE DATA STRUCTURE
	// DO NOT MODIFY IT
	// SUPPRESSING THE SCM CLASS OBJECT WILL DESTROY THE RESULTS	
	std::vector<std::vector<int>>* returnAttractors() {
		return &periodicCells;
	}
	std::vector<std::vector<int>>* returnDomainsOfAttraction() {
		return &domainsOfAttraction;
	}


protected:
	// INHERIT FUNCTIONS FROM PARENT CLASS
	using cellMap::cellMap;
	using cellMap::getCellImage; 
	using cellMap::getCellImageRefined; 
	using cellMap::nbCells;
	using cellMap::sinkCell;


	// COMPUTE IMAGE OF ALL CELLS OF cellMap
	void computeMapImage() {
		cellMapImage[sinkCell] = sinkCell;
		#pragma omp parallel for default(shared) if(UseParallel!=0) 
		for (int it=0; it<nbCells-1; it++) {
			double time=0.0;
			double end = 0.0;
			double begin = 0.0;
			auto begin1=std::chrono::steady_clock::now();	
			cellMapImage[it] = getCellImage(it);
			auto end1=std::chrono::steady_clock::now();
			auto timei=std::chrono::duration_cast<std::chrono::nanoseconds>(end1-begin1).count();
			time = timei;
			/*
			if (it%10000==0) {
			#if _OPENMP
			#pragma omp critical
				{
					cout<<"end and begin: "<<end<<" "<<begin<<endl;
					cout<<"cell number: "<<it<<" computed by thread: "<<omp_get_thread_num()<<" in "<<time<<" nanoseconds"<<endl;
				}
			#else
				cout<<"cell number: "<<it<<" computed in "<<timei<<" milliseconds"<<endl;
			#endif
			}*/
		}
	}	

	
	//TELLS IF CELL cCell HAS AN ATTRACTOR (YET) OR NOT
	bool hasAttractor(const int cCell) {	
		for (int it=0; it<attractors.size(); it++) {
			if (cCell==attractors[it])
				return true;	
		}
		return false;
	}


	// EXPLORES cellMapImage USING TARJAN'S ALGORITHM TO FIND ALL ATTRACTORS AND THEIR DOMAIN OF ATTRACTION
	void findAttractors() {
		int it3, cell; // Number of periodic group

		std::vector<int> currentPoints;
		for (int it=0; it<nbCells-1; it++) {
			cell=it;
			if (hasAttractor(cellMapImage[cell])==false) {
				//keep exploring untill finding a circularity (ie new attractor) or finding an already existing attractor
				while (cellMapImage[cell]!=nbCells+2 
				       && hasAttractor(cellMapImage[cell])==false) {
					currentPoints.push_back(cell);
					cell=cellMapImage[cell];
					cellMapImage[currentPoints[currentPoints.size()-1]]=nbCells+2;
				}
				// Define in which attractor the sequence of points is
				if (cellMapImage[cell]==nbCells+2) { //new attractor found
					attractors.push_back(cell);  
					for (size_t it2=0; it2<currentPoints.size(); it2++)
						cellMapImage[currentPoints[it2]]=cell;

					if (cell!=currentPoints[currentPoints.size()-1]) { // if there is a cycle
						// Find cells belonging to the cycle within currentPoints
						// Namely elements ranging from it3 to currentPoints.size() in currentPoints 
						it3=currentPoints.size()-1;
						while(cell!=currentPoints[it3])
							it3-=1;
						std::vector<int>::const_iterator first = currentPoints.begin() + it3;
						std::vector<int>::const_iterator last = currentPoints.end();
						std::vector<int> attractor(first, last);
						periodicCells.push_back(attractor);
					}
					else { // if there is a stable equilibrium point
						std::vector<int> attractor(1);
						attractor[0] = cell;
						periodicCells.push_back(attractor);
					} 
				} else { //already existing attractor found
				for (size_t it2=0; it2<currentPoints.size(); it2++)
					cellMapImage[currentPoints[it2]]=cellMapImage[cell];
				}
		}
			currentPoints.clear();
		}
	}


	void findDomainsOfAttraction() {
		bool found;
		int it2;
		domainsOfAttraction.resize(attractors.size());
		for (int it=0; it<nbCells; it++) {
			it2=1;
			found=false;
			while(it2<attractors.size() && found==false) {
				if (cellMapImage[it]==attractors[it2] ) {
					domainsOfAttraction[it2].push_back(it);
					found=true;
				}
				it2+=1;
			}
		}
	}

	// TESTS IF CELL cell IS IN ATTRACTOR attractorId
	bool isInAttractor(int cell, int attractorId) {
		for (int it=0; it<periodicCells[attractorId].size(); it++) {
			if (cell==periodicCells[attractorId][it]) {
				return true;
			}
		}
		return false;
	}
	
	

	// MERGE TWO ATTRACTORS
	// Could be made more readable by using iterators instead of int for attractor1 and attractor2....
	void mergeAttractors(int attractor1, int attractor2) {
		
		// Merge attractors
		periodicCells[attractor1].insert(periodicCells[attractor1].end(), periodicCells[attractor2].begin(), periodicCells[attractor2].end());
		
		// Merge Domains Of Attraction
		if (computeDomainOfAttraction==true) {
			domainsOfAttraction[attractor1].insert(domainsOfAttraction[attractor1].end(), domainsOfAttraction[attractor2].begin(), domainsOfAttraction[attractor2].end());
		}
	}
	int* cellMapImage;				// Image of all cells of the cellMap
	bool computeDomainOfAttraction;			// Compute domains of attraction or only attractors
	std::vector<int> attractors;			// ID of all attractors
	std::vector<std::vector<int>> periodicCells;	// Members of all attractors	
	std::vector<std::vector<int>> domainsOfAttraction;		// Domain of attraction of all attractors

};

#endif
