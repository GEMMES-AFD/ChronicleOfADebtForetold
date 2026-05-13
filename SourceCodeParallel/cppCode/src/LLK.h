#ifndef MINDIST_h
#define MINDIST_h

#if UseRCPP==1
	#include "preprocRCPP_R.h"
#else
	#include "preproc.h"	
#endif

#include "ODE.h"

using namespace std;


// TO DO: pass yInit and parmslength in the model definition to allow for more abstraction
// TO DO: specific class for varMinDist to improve code readability
// TO DO: use other norm than absolute value (norm_2 or norm_inf) => main priority

// COMPUTES DISTANCE (objective function for calibration) 
template<typename T>
class LLK { 
public:
	// CONSTRUCTOR
	LLK(ODE<T>* aModelIn, const int parmsLengthIn,
		T* yInitIn, int nVarMinDistIn, int* varMinDistIn, T* dataMinDistIn, 
		T* samplingTimeMinDistIn, T* pointsWeightMinDistIn, 
		int* nObsVarMinDistIn) :
		myModel(aModelIn), parmsLength(parmsLengthIn),
		yInit(yInitIn), nVarMinDist(nVarMinDistIn), varMinDist(varMinDistIn), dataMinDist(dataMinDistIn),
	        samplingTimeMinDist(samplingTimeMinDistIn), pointsWeightMinDist(pointsWeightMinDistIn),
		nObsVarMinDist(nObsVarMinDistIn) {}

	T Evaluate(T* parms) {
		T out = 0;
		try {
			out = EvaluateWithoutTryCatch(parms);
		}
		catch(...) {
			#if VerboseCMAES>2
				#pragma omp critical
					{
						#if UseRCPP==1
							Rcpp::Rcout<<"WARNING: Failure to estimate objective function: unknown error occurred. For parameters: ";
							for (int it=0; it<parmsLength; it++) {
								Rcpp::Rcout<<parms[it]<<" ";
							}
							Rcpp::Rcout<<". Setting objective function to 1e50 and continuing."<<endl;
						#else
							cout<<"WARNING: Failure to estimate objective function: unknown error occurred. For parameters: ";
							for (int it=0; it<parmsLength; it++) {
								cout<<parms[it]<<" ";
							}
							cout<<". Setting objective function to 1e50 and continuing."<<endl;
						#endif
				}
			#endif
					out = 1e50;
		}
		if (std::isnan(out)) {
			#if VerboseCMAES>2
				#pragma omp critical
					{
						#if UseRCPP==1
							Rcpp::Rcout<<"WARNING: Failure to estimate objective function: value returned is NaN. For parameters: ";
							for (int it=0; it<parmsLength; it++) {
								Rcpp::Rcout<<parms[it]<<" ";
							}
							Rcpp::Rcout<<". Setting objective function to 1e50 and continuing."<<endl;
						#else
							cout<<"WARNING: Failure to estimate objective function: value returned is NaN. For parameters: ";
							for (int it=0; it<parmsLength; it++) {
								cout<<parms[it]<<" ";
							}
							cout<<". Setting objective function to 1e50 and continuing."<<endl;
						#endif
				}
			#endif
			out = 1e50;
		}
		return out;
	}

protected:
	// UPDATE FIT AT ITERATION it IF NEEDED 
	T updateFit(const int it, int* currentPosVarMinDist, T* trajectory) {
		T out = 0;
		for (int it1=0; it1<nVarMinDist; it1++) {
			if (abs(myModel->getTInit() + it*myModel->getHOut()- samplingTimeMinDist[currentPosVarMinDist[it1]])<myModel->getHOut()/2) {
				if (dataMinDist[currentPosVarMinDist[it1]]!=0.0) {
					out+=pointsWeightMinDist[currentPosVarMinDist[it1]]*std::abs((getSimulatedValue(trajectory, it, it1) - dataMinDist[currentPosVarMinDist[it1]])/dataMinDist[currentPosVarMinDist[it1]]); //Norm 1
					//out+=pointsWeightMinDist[currentPosVarMinDist[it1]]*pow((getSimulatedValue(trajectory, it, it1) - dataMinDist[currentPosVarMinDist[it1]])/dataMinDist[currentPosVarMinDist[it1]], 2); //Norm 2
				} else {
					out+=pointsWeightMinDist[currentPosVarMinDist[it1]]*std::abs((getSimulatedValue(trajectory, it, it1) - dataMinDist[currentPosVarMinDist[it1]]));
				}
			currentPosVarMinDist[it1]++;
			}
		}
		return out;
	}

	T getSimulatedValue(const T* trajectory, const int it, const int it1) {
		#if ReturnRK4==2
			return trajectory[myModel->getNRowOut()*it + varMinDist[it1]];
		#elif ReturnRK4==3
			if(varMinDist[it1]<myModel->getNV()) {
				return trajectory[myModel->getNRowOut()*it + varMinDist[it1]];
			} 
			else {
				return trajectory[myModel->getNRowOut()*it + varMinDist[it1]+myModel->getNV()];
			}
		#else
			throw std::exception("ReturnRK4 must be set to 2 or 3 to compute LLK !");
		#endif

	}
	// COMPUTES DISTANCE FOR A SET OF PARAMETERS VALUES (no error management)	
	T EvaluateWithoutTryCatch(T* parms) {
		// INIT
		T fit = 0.0;
		T* trajectory = new T[myModel->getNRowOutSolve()];
		int currentPosVarMinDist[nVarMinDist];
		// To store number of points already explored in dataMinDist
		currentPosVarMinDist[0] = 0;
		for (int it1=1; it1<nVarMinDist; it1++) {
			currentPosVarMinDist[it1] = currentPosVarMinDist[it1-1]+nObsVarMinDist[it1-1];
		}	
		// COMPUTE TRAJECTORY
		myModel->solve(yInit, parms, trajectory);
		for (int it=0; it<myModel->getNt(); it++) {
			fit+=updateFit(it, currentPosVarMinDist, trajectory);
		}
		free(trajectory);
		return fit;
	}

	ODE<T>* myModel;
	const int parmsLength;
	T* yInit;
	int nVarMinDist;
	int* varMinDist;
	int* nObsVarMinDist;
	T* dataMinDist;
	T* samplingTimeMinDist;
	T* pointsWeightMinDist;
//	T* xyTrajectory;
};


#endif
