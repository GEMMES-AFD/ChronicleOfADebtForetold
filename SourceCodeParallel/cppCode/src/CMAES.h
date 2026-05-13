#ifndef CMAES_H
#define CMAES_H

#include <armadillo>

#if UseRCPP==1
	#include "preprocRCPP_R.h"
#else
	#include "preproc.h"	
#endif

using namespace arma; 

/******************************/
/*  Optimization using CMAES  */
/******************************/
// evaluateClass is a class defining the objective function for optimization, it must have the following structure:
// class evaluateClass {
// public:
// 	T Evaluate(arma::Col<T>);          // The function to minimize
// };


// TO DO: cut main function (optimize) in multiple functions to make code structure explicit
// TO DO: find way to have one single parallel block to save computation time
// TO DO: if point above is impossible, find way to set random seed only once to save computation time

template<typename T, class evaluateClass>
class CMAES {
public:
	// CONSTRUCTOR
	CMAES(evaluateClass& anEvaluateClass, 
	      int lambdaIn,
	      double sigmaIn, 
	      int nIterMaxIn, 
	      double tolIn,
	      T* mLowerIn=nullptr, 
	      T* mUpperIn=nullptr, 
		  int nReDrawMaxIn=100,
	      bool useParallelIn=true) : 
	        myEvaluateClass(anEvaluateClass),
	        lambda(lambdaIn), sigma(sigmaIn), nIterMax(nIterMaxIn), 
       	       	tol(tolIn), useParallel(useParallelIn), normalize(mLowerIn!=nullptr),
       		mLower(mLowerIn), mUpper(mUpperIn), nReDrawMax(nReDrawMaxIn) {};



	// WRAPPER AROUND evaluateClass::Evaluate
	// TO NORMALIZE PARAMETERS (IF NEEDED) AND PENALIZE IF OUT OF BOUNDS
	T EvaluateBounded(arma::Col<T> m) {
		T out;
		if (normalize) {
			unnormalizeM(m);
			// add penalty
			bool outOfBounds = false;
			for (int it=0; it<m.n_elem; it++) {
				if (m[it]<mLower[it] || m[it]>mUpper[it]) {
					#if VerboseCMAES>0
						#pragma omp critical
							{
							Rcpp::Rcout<<mLower[it]<<" "<<m[it]<<" "<<mUpper[it]<<endl;
								#if UseRCPP==1
										Rcpp::Rcout<<"WARNING: Parameter number "<<it<<" out of exploration domain. Setting objective function to 1e50 and continuing."<<endl;
								#else
										cout<<"WARNING: Parameter number "<<it<<" out of exploration domain. Setting objective function to 1e50 and continuing."<<endl;
								#endif
							}
					#endif
					outOfBounds = true;
				}
			}
			if (outOfBounds==true) {
				out=1e50; // Add arbitrarily large value as penalty
			}
			else {
				out=myEvaluateClass.Evaluate(m);
			}
		} 
		else {
			out=myEvaluateClass.Evaluate(m);
		}
		return out;
	}

	// SEARCH THE VALUES OF m THAT MINIMIZE evaluateClass.Evaluate(.) USING CMA-ES ALGORITHM
	double Optimize(arma::Col<T>& m) {
		
		// INITIALIZATION //
		if(normalize) normalizeM(m);

		// INITIALIZING PARMS FOR CMAES. SEE HANSEN(2016) FOR MORE INFORMATION
		const int n = m.n_elem;
		double mu = floor(lambda/2);
		arma::Col<T> w = createw(lambda, n);
		double muEff = pow(accu(w(span(0, mu-1))), 2)/accu(pow(w(span(0, mu-1)), 2)); 
		double EN01 = std::pow(n, 0.5)*(1-1/(4*n)+1/(21*std::pow(n, 2)));
		double csig = (muEff+2)/(n + muEff+5);
		double dsig = 1 + 2*std::max(0.0, std::pow((muEff-1)/(n+1), 0.5)-1) + csig;
		double cc = (4 + muEff/n)/(n + 4 + 2*muEff/n);
		double c1 = 2/(std::pow(n + 1.3, 2) + muEff);
		double cb = 1;
		double cmu = std::min(1-c1, 2*(muEff- 2 + 1/muEff)/(std::pow(2+n, 2) + 2*muEff/2));
		double hsig = 0;
		int indexPastVal = 0;
		int maxPastVal = 95 + floor(30*n/lambda);
		int indexNewestVal = 0;
		arma::Col<T> NewestVal(25);
		arma::Col<T> OldestVal(25);

		int nTrial;
		arma::Col<T> yw(n);
		arma::Mat<T> C = eye(n, n);
		arma::Mat<T> zlambda(n, lambda);
		arma::Mat<T> ylambda(n, lambda);
		arma::Mat<T> xlambda(n, lambda);
		arma::Col<T> flambda(lambda);
		arma::Col<T> psig = zeros(n);
		arma::Col<T> pc = zeros(n);
		arma::uvec flambdaIndex(lambda);
		arma::Col<T> D(n);
		arma::Mat<T> B(n, n);
		arma::Col<T> tempM;
		arma::Col<T> pastVal(maxPastVal);
		arma::Col<int> pastAvgResults(n, fill::zeros);
		arma::Col<T> wo = w;
		int indexPastAvgResults = 0;
		/*************/
		// MAIN LOOP //
		/*************/
		
		for (int i=0; i<nIterMax; ++i) {
#if VerboseCMAES>0
	#if UseRCPP==1
					Rcpp::Rcout<<"Beginning iteration number "<<i<<endl;
	#else
					cout<<"Beginning iteration number "<<i<<endl;
	#endif
#endif
#if VerboseCMAES>2
		arma::Col<T> mForPrint(n);
		for (int itVerbose=0; itVerbose<n; itVerbose++) mForPrint(itVerbose) = m(itVerbose);
		if (normalize) unnormalizeM(mForPrint);
		#if UseRCPP==1
						Rcpp::Rcout<<"Value of mean point is: \n"<<mForPrint.t()<<endl;
		#else
						cout<<"Value of mean point is:  \n"<<mForPrint.t()<<endl;
		#endif
#endif
			// CREATE NEW GENERATION OF CANDIDATES AND EVALUATE THEIR FITNESS
			eig_sym(D, B, C);
			D = pow(D, 0.5);
			// Note: we can't have a parallel construct nested in a single construct,
			//       so that we need to start a new parallel region at each iteration of the for loop

			#pragma omp parallel default(shared) private(nTrial, tempM) if(useParallel)
		       	{
				arma_rng::set_seed_random();
				#pragma omp for
				for (int j=0; j<lambda; ++j) { //Generate offsprings and compute their fitness
					nTrial = 0;
					do {
						zlambda.col(j) = randn(n);
						ylambda.col(j) = B*diagmat(D)*zlambda.col(j);
						xlambda.col(j) = m + sigma*ylambda.col(j);
						tempM = vec(xlambda.colptr(j), n, false, false);
						flambda(j) = EvaluateBounded(tempM);
						nTrial++;
					} while (flambda(j)>=1e49 && nTrial<nReDrawMax); // Re-generate points until the LLK is successfully computed
					#if VerboseCMAES>0
						if (nTrial==nReDrawMax) {
							#pragma omp critical
								{
									#if UseRCPP==1
											Rcpp::Rcout<<"WARNING: Maximum number of redraw reached, for point "<<j<<". Setting objective function to 1e50 and continuing."<<endl;
									#else
											cout<<"WARNING: Maximum number of redraw reached, for point "<<j<<". Setting objective function to 1e50 and continuing."<<endl;
									#endif
							}
						#if VerboseCMAES>1
						} else if (nTrial>1) {
							#pragma omp critical
								{
									#if UseRCPP==1
											Rcpp::Rcout<<"WARNING: Point succesfully estimated after "<<nTrial-1<<" redraws."<<endl;
									#else
											cout<<"WARNING: Point succesfully estimated after "<<nTrial-1<<" redraws."<<endl;
									#endif
								}
						#endif
						}
					#endif
				}
			}
			flambdaIndex = sort_index(flambda, "ascend"); // sort points by fitness
			yw = ylambda.cols(flambdaIndex(arma::span(0, mu-1)))*w(span(0, mu-1));
			//cout<<flambdaIndex(0)<<" "<<flambda(flambdaIndex(0))<<endl;
			if (flambda(flambdaIndex(0))>=1e50) {
				#if UseRCPP==1
					Rcpp::Rcout<<"Failure: Calculation of the objective function failed (i.e. returned value >= 1e50) for all points explored, at iteration "<<i<<".Early code termination."<<std::endl;
				#else
					cout<<"Failure: Calculation of the objective function failed (i.e. returned value >= 1e50) for all points explored, at iteration "<<i<<". Early code termination."<<std::endl;
				#endif
				break;
			}

			if (arma::min(flambda(flambdaIndex(arma::span(0, mu-1))))>=1e50) {
				#if UseRCPP==1
					Rcpp::Rcout<<"WARNING: Calculation of the objective function failed (i.e. returned value >= 1e50) for more than 50% of points explored, at iteration "<<i<<" this might prevent convergence from occurring. Continuing nonetheless."<<std::endl;
				#else
					cout<<"WARNING: Calculation of the objective function failed (i.e. returned value >= 1e50) for more than 50% of points explored, at iteration "<<i<<" this might prevent convergence from occurring. Continuing nonetheless."<<std::endl;
				#endif
			}

			// SELECTION AND RECOMBINATION
			m = m + cb*sigma*yw; //update m taking the mean of the selected points
			// STEP-SIZE CONTROL
			psig = (1 - csig)*psig + std::pow(csig*(2-csig)*muEff, 0.5)*B*zlambda.cols(flambdaIndex(span(0, mu-1)))*w(span(0, mu-1));
			sigma = sigma*exp(csig/dsig*(norm(psig)/EN01-1));

			//COVARIANCE MATRIX ADAPTATION
			for (int j=mu;j<lambda; ++j) {
				wo(j)= w(j)*n/pow(norm(B*zlambda.col(flambdaIndex(j)), 2), 2);
			}
			hsig = ((norm(psig)/std::pow(1 - std::pow(1.0-csig, 2.0*(n+1.0)), 0.5)) < ((1.4 + 2/(n+1))*EN01)) ? 1.0 : 0.0;
			pc = (1-cc)*pc + hsig*std::pow(cc*(2-cc)*muEff, 0.5)*ylambda.cols(flambdaIndex(arma::span(0, mu-1)))*w(span(0, mu-1));
			C = (1 + c1*(1-hsig)*cc*(2-cc) - c1 - cmu*accu(wo))*C + c1*pc*pc.t() + cmu*ylambda.cols(flambdaIndex)*diagmat(wo)*ylambda.cols(flambdaIndex).t();
			C = trimatu(C) + trimatl(C, -1); //Enforce symmetry of C
			
			
			//if (std::abs((flambda(flambdaIndex(0))-flambda(flambdaIndex(mu-1)))/flambda(flambdaIndex(0)))<tolFitness) { //check for flatness
			//  sigma = sigma*exp(0.2+csig/dsig);
			//  Rcpp::Rcout<<"Warning: Flat fitness function"<<std::endl;
			//}
				
			// TERMINATION CRITERIA
			// There are multiple stopping criterias,
			// see Hansen(2016) for more information
			if (max(pc)<tol && max(sigma*C.diag())<tol) {
				#if UseRCPP==1 // Different function for print depending on if calling CMAES from c++ or from R
				 Rcpp::Rcout<<"Terminating successfully: pc and C's values are below threshold, at iteration: "<<i<<std::endl;
				#else
				 cout<<"Terminating successfully: pc and C's values are below threshold, at iteration: "<<i<<std::endl;
				#endif
				break;
			}
			if (max(D)>pow(10, 14)*min(D)) {
				#if UseRCPP==1
				 Rcpp::Rcout<<"Terminating successfully: Excessive condition number of the covariance matrix, at iteration: "<<i<<std::endl;
				#else
				 cout<<"Terminating successfully: Excessive condition number of the covariance matrix, at iteration: "<<i<<std::endl;
				#endif
				break;
			}

			(indexPastVal<maxPastVal-1) ? indexPastVal++ : indexPastVal = 0;
			(indexNewestVal<25-1) ? indexNewestVal++ : indexNewestVal=0;
			OldestVal(indexNewestVal) = pastVal(indexPastVal);
			pastVal(indexPastVal) = EvaluateBounded(m);
			NewestVal(indexNewestVal) = pastVal(indexPastVal);
			if (i>maxPastVal+25 && std::abs(mean(OldestVal) - mean(NewestVal))<tol & std::abs(median(OldestVal) - median(NewestVal))<tol) {
				#if UseRCPP==1
				 Rcpp::Rcout<<"Terminating successfully: Stagnation of function value at m, at iteration: "<<i<<std::endl;
				#else
				 cout<<"Terminating successfully: Stagnation of function value at m, at iteration: "<<i<<std::endl;
				#endif
				break;
			}

			(indexPastAvgResults<n-1) ? indexPastAvgResults++ : indexPastAvgResults = 0;
			flambda(flambdaIndex(mu-1)) - flambda(flambdaIndex(0))<tol ? pastAvgResults(indexPastAvgResults)=1 : pastAvgResults(indexPastAvgResults) = 0;
			if (accu(pastAvgResults)>n/3) {
				#if UseRCPP==1
				 Rcpp::Rcout<<"Terminating successfully: Constant value for all best mu points, at iteration: "<<i<<std::endl;
				#else
				 cout<<"Terminating successfully: Constant value for all best mu points, at iteration: "<<i<<std::endl;
				#endif
				break;
			}
			/*	
			if (norm(0.1*sigma*D(n-1)*B.col(n-1)/m)<tol) {
				cout<<"tol: "<<tol<<" condition: "<<norm(0.1*sigma*D(n-1)*B.col(n-1)/m)<<endl;
				#if UseRCPP==1
				 Rcpp::Rcout<<"Terminating successfully: A shock on the principal axis of C does not significantly change m, at iteration: "<<i<<std::endl;
				#else
				cout<<"Terminating successfully: A shock on the principal axis of C does not significantly change m, at iteration: "<<i<<std::endl;
				#endif
				break;
			}
			*/
			//if (flambda(flambdaIndex(mu-1))>1e45) break;
			if (i==nIterMax-1) {
			#if UseRCPP==1
				Rcpp::Rcout<<"Terminating before convergence: maximum number of iteration reached. Returning last explored point."<<std::endl;
			#else
				cout<<"Terminating before convergence: maximum number of iteration reached. Returning last explored point."<<std::endl;
			#endif
			}
		}
		// generate and return output (best solution m and value at the best solution)
		if(normalize) unnormalizeM(m);
			
		return pastVal(indexPastVal);
	}
	T Evaluate(arma::Col<T> m) {
		return myEvaluateClass.Evaluate(m);
	}
protected:
	evaluateClass& myEvaluateClass;
	int lambda; 			// number of candidates per generation
	double sigma;			// initial overall dispersion of candidates
	int nIterMax;			// maximum number of iterations
	double tol;			// tolerance for convergence criterion
	double tolFitness = 1e-3;	// tolerance for convergence criterion
	bool useParallel;		// run algorithm in parallel or not
	T* mLower;
	T* mUpper;
	bool normalize;
	int nReDrawMax;

	// SET ALL PARMS TO VALUES BETWEEN 0 AND 1 DEPENDING ON THEIR LOWER AND UPPER BOUNDS (mLower AND mUpper) 
	void normalizeM(arma::Col<T>& m) {
		for(int it=0; it<m.n_elem; it++) {
			m[it]-=mLower[it];
			m[it]/=(mUpper[it]-mLower[it]);
		}
	}
	
	// RECOVER REAL PARAMETERS VALUE FROM THEIR STANDARDIZED VALUE
	void unnormalizeM(arma::Col<T>& m) {
		for (int it=0; it<m.n_elem; it++) {
			m[it]*=(mUpper[it]-mLower[it]);
			m[it]+=mLower[it];
		}
	}
	
	// INITIALIZE VECTOR w
	arma::Col<T> createw(const int LAMBDA, const int N) { 
		const int MU = floor(LAMBDA/2);
		arma::Col<T> out(LAMBDA);
		for (int i=0; i<LAMBDA; ++i) {
			out(i) = log((LAMBDA+1)/2) - log(i+1);
		}
		double MUEFF = pow(accu(out(span(0, MU-1))), 2)/accu(pow(out(span(0, MU-1)), 2));
		double MUEFFminus = pow(accu(out(span(MU, LAMBDA-1))), 2)/accu(pow(out(span(MU, LAMBDA - 1)), 2));
		double C1 = 2/(std::pow(N + 1.3, 2) + MUEFF);
		double CMU = std::min(1-C1, 2*(MUEFF- 2 + 1/MUEFF)/(std::pow(2+N, 2) + 2*MUEFF/2));
		double ALPHAMUminus = 1 + C1/CMU;
		double ALPHAMUEFFminus = 1 + 2*MUEFFminus/(MUEFF+2);
		double ALPHAPOSDEFminus = (1 - C1 - CMU)/(N*CMU);
		double SUMTEMPplus = accu(out(span(0, MU-1)));
		double SUMTEMPminus = std::abs(accu(out(span(MU, LAMBDA-1))));
		for (int i = 0; i<MU;++i) {
			out(i) *= 1/SUMTEMPplus;
		}
		for (int i = MU; i<LAMBDA;++i) {
			out(i) *= std::min(std::min(ALPHAMUminus, ALPHAMUEFFminus), ALPHAPOSDEFminus)/SUMTEMPminus;
		}
		return out;
	}

	
};

#endif
