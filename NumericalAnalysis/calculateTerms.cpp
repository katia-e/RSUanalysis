#include <iostream> // Perturbation using Jacobi iterations
#include <stdio.h>
#include <stdlib.h>
#include "stdint.h"
#include <sys/time.h>
#include <cmath>
using namespace std;
#include "libStateSpace.h"
//#include <iomanip>
#include <sstream>


string toStringPrecision(double x, int precision){
	ostringstream out;
	out<<std::setprecision(precision)<<std::fixed<<x;
	return out.str();
}
string toStringScient(double x, int precision){
	ostringstream out;
	out<<std::setprecision(precision)<<std::scientific<<x;
	return out.str();
}

int main(int argc,char *argv[])
{
	/* Default parameters */
	
	/* System parameters*/
	int K = 7;					// max number of users in the system
	int C = 3;					// buffer capacity
	int M = 3;					// number of ranges	
	double beta[] = {1.0/30, 1.0/15, 1.0/30}; //  transition rate between ranges
	double data_rate[] = {0.1, 1.0, 0.1};
	double car_rate = 0.1;		// rate of car arrivals	
	double lam = 1;				// packet arrival rate	
	
	/* Analysis parameters */
	int nTerms = 10;					// Number of terms
	int nIter = 100;					// Number of iterations for each term
	double stopping = 0.0001; 			// threshold of the spotting creteria 	
	double MuJac = 0;					// reference point for series expansion	
	double W = 0.9; 					// Stabilizing parameter
	string fOutput = "";		// output file for the result
	
	/* Parse comand line parameters if given */
	for(int i = 1; i<argc; i++) 
		switch (i){			
			case 1: MuJac = atof(argv[i]); 
			case 2: car_rate = atof(argv[i]);
			case 3: K = atoi(argv[i]); 
			case 4: C = atoi(argv[i]); 
			case 5: M = atoi(argv[i]); 
			case 6: lam = atof(argv[i]);
			case 7: nTerms = atoi(argv[i]);
			case 8: stopping = atof(argv[i]);
			case 9: W = atof(argv[i]);
			case 10: beta[0] = atof(argv[i]);
			case 11: beta[1] = atof(argv[i]);
			case 12: beta[2] = atof(argv[i]);
			case 13: data_rate[0] = atof(argv[i]);
			case 14: data_rate[1] = atof(argv[i]);
			case 15: data_rate[2] = atof(argv[i]);	
			case 16: fOutput = argv[i];
		};
	double percision = 1;
	if (fOutput.length() == 0)
		/* Generate distinct name of the output file */	
		fOutput = "QTermsK"+to_string(K)+
						"-M"+to_string(M)+"-C"+to_string(C)+
						"-mu"+toStringPrecision(MuJac,percision)+"-lam"+toStringPrecision(lam,percision)+
						"-W"+toStringPrecision(W,percision) + "-stop"+toStringScient(stopping,0)+".dat";
	cout<<endl<<"Mu0 = "<<MuJac<<endl
				<<"car_rate = "<<car_rate<<endl
				<<"K = "<<K<<endl
				<<"C = "<<C<<endl
				<<"M = "<<M<<endl
				<<"Beta = "<<beta[0]<<", "<<beta[1]<<", "<<beta[2]<<endl
				<<"lambda = "<<lam<<endl
				<<"nTerms = "<<nTerms<<endl
				<<"eps = "<<stopping<<endl
				<<"W = "<<W<<endl
				<<"Output file is "<<fOutput<<endl
				<<"beta = "<<beta[0]<<beta[1]<<beta[2]<<endl
				<<"data_rate = "<<data_rate[0]<<data_rate[1]<<data_rate[2];

	/* Parameters of the state space */
	int DIM = M*(C+1);				// Number of dimentions in the queueing model
	cout<<"\nDimentions "<<DIM;	
	
	int nStates = binomial(K, DIM); // Number of states in the queueing model
	cout<<"\nState space size "<<nStates;

	/* StateSpace defines structure for the state space. 
		Dimenssion of the array @StateSpace is defined by @DIM and has to be hardcoded.
			It is neceessary to check that @DIM equals the depth of @StateSpace.
			Function allocateStateSpace() supports values of DIM = 6, 10, 12.
		When changing @M or @C, check the value of DIM, 
			uncomment appropriate initialisation for StateSpace and recompile.
	*/
//	int**********StateSpace;	StateSpace = new int********* [K+1];	// DIM = 10
//	int******StateSpace;	StateSpace = new int***** [K+1]; // DIM = 6
	int************StateSpace;	
	StateSpace = new int*********** [K+1]; // DIM = 12
	allocateStateSpace(StateSpace, K);	
	
	int iState = 0;
	int size = 0;

	
	// define vector of data rates mu
	vect mu;
	mu.init(DIM-M);
	mu.fill(data_rate, M);
	cout<<"\n Data rates in zones: ";
	print(mu);	
	
	//getchar();	
	// Predefinitions for calculation
	vect *ssdTermsJac; ssdTermsJac = new vect [nTerms];
	vect QTerms; QTerms.init(nTerms);
	vect nCarsTerms; nCarsTerms.init(nTerms);
	for(int i=0; i<nTerms; i++)
		ssdTermsJac[i].init(nStates);
	stState state, jstate; 
	state.init(DIM); jstate.init(DIM);
	timeval dw1, dw2;
	for(int iTerm = 0; iTerm < nTerms; iTerm++){// start calculating terms
		gettimeofday(&dw1, 0);
		bool it_is_not_first_term = (iTerm>0);
		cout<<"\n"<<iTerm<<"th term";
		bool CONTINUE = true;
		int iIter = 0;
		double QTermPrev = 0;
		while(CONTINUE){
			iIter++;			
			QTerms.x[iTerm] = 0;
//		for(int iIter = 1; iIter<nIter; iIter++){// start Jacobi iterations
			//bool it_is_last_iteration = (iIter==nIter-1);
			state.to0(); state.st[0] = K;
			for(int iState = nStates-1; iState>0; iState--){ //calculation for each state
				int nCars = state.sum();
				double Pr = 0; double bi = 0; double aii = 0;
				double tmp = 0; int x; int from, to,up, down; // sup variable
				int jState;
				//  All possible transissions from iState a11
				// Get from iState by a car arrival
				if(nCars < K){ // then car arrival is possible
					aii -= car_rate;
				}
				// Transitions with packet arrivals and car arrivals
				for(int n = 0; n < C+1; n++)					
					for(int m = 0; m < M; m++){
						x = state.st[M*n+m];
						if(n<C)  
							aii -= x*(lam+beta[m]);
						else           
							aii -= x*beta[m];
					}
				//  Get from iState by a packet transmission
				if(state.sum(M,DIM)>0){
					tmp = 0; 
					for(int n=1;n<C+1;n++)
						for(int m=0;m<M;m++){
							int ind = M*n+m;
							tmp += state.st[ind]*mu.x[ind-M];
						}
					aii -= MuJac*tmp;
                    if(it_is_not_first_term)// finding states where system hops by packet transmission
                         bi += ssdTermsJac[iTerm-1].x[iState]*tmp;                          
				}
				// Get to iState by a car arrival
				if(state.st[0]>0){ // then getting to iState by a car arrival is possible
					jstate.copy(state);
					jstate.st[0] --; // % From this state	
					jState = getStateNum(StateSpace, jstate);				
					Pr += car_rate*ssdTermsJac[iTerm].x[jState];
				}
				// Get to iState by a packet arrival   
				if(state.sum(M,DIM)>0){// then getting to iState by packet arrival is possible
					tmp = 0;
					for(int n=1;n<C+1;n++){					
						for(int m=0;m<M;m++){
							up = M*n+m;
							if(state.st[up]>0){
								down = M*(n-1)+m;
								x = state.st[down]+1;
								jstate.copy(state);
								jstate.st[up] --; jstate.st[down] ++; 
								jState = getStateNum(StateSpace, jstate);	
								tmp+= x*ssdTermsJac[iTerm].x[jState];
							}
						}
					}
					Pr += lam*tmp;
				}
				// Get to iState by switching the range
				tmp = 0;
				for(int n=0;n<C+1;n++){					
					for(int m=1;m<M;m++){
						up = M*n+m;
						if(state.st[up]>0){
							down = M*n+m-1;
							x = state.st[down]+1;
							jstate.copy(state);
							jstate.st[down]++; jstate.st[up]--;																			
							jState = getStateNum(StateSpace, jstate);							
							tmp += x*ssdTermsJac[iTerm].x[jState]*beta[m-1];
						}
					}
				}
				Pr += tmp;          
				// Get to iState by one of the cars leaving the system
				if(state.sum()<K) {//then reaching this state is possible by a car departure
					tmp = 0;
					for(int n=1;n<C+2;n++){					
						up = M*n-1;
						x = state.st[up]+1;
						jstate.copy(state);
						jstate.st[up] ++;					
						jState = getStateNum(StateSpace, jstate);	
						tmp += x*ssdTermsJac[iTerm].x[jState];
					}
					Pr += beta[M-1]*tmp;
				}
				  // Get to iState by transmitting a packet
				if(state.sum(0,DIM-M)>0){ // then getting to the currect state by packet arrival is possible
					 tmp = 0;
					 for(int n = 1; n < C+1; n++){					
						for(int m = 0; m < M; m++){
							down = M*(n-1)+m; 
							if(state.st[down]>0){
								up = M*n+m;
								x = state.st[up]+1;
								jstate.copy(state);
								jstate.st[down]--; jstate.st[up]++;
								//jState = getStateNum12dim(StateSpace, jstate);
								jState = getStateNum(StateSpace, jstate);
								double mu_scheduled = mu.x[up-M];
								tmp += x*MuJac*mu_scheduled*ssdTermsJac[iTerm].x[jState];//flag
                                if(it_is_not_first_term) // Perturbed part
                                    bi -= x*mu_scheduled*ssdTermsJac[iTerm-1].x[jState]; 
							}
						}
					 }
					Pr = Pr + tmp; 
				}			       
//				ssdTmp.x[iState] = (bi-Pr)/aii;; // Store previous value
//				ssdTermsJac[iTerm].x[iState] = (bi-Pr)/aii;  
				ssdTermsJac[iTerm].x[iState] =  (1-W)*ssdTermsJac[iTerm].x[iState]+W*(bi-Pr)/aii;				
				int nPacketsInStates = 0;				
				for(int n = 1; n<C+1; n++)
					nPacketsInStates += state.sum(n*M,(n+1)*M)*n;
				QTerms.x[iTerm] += nPacketsInStates*ssdTermsJac[iTerm].x[iState];	
				
//				nCarsTerms.x[iTerm] += state.sum()*Pr;
				dcrvect(state, state, K);
			} // all states are calculated
			// normalize terms
            if(it_is_not_first_term) // sum(term) = 0
                ssdTermsJac[iTerm].x[0] -= ssdTermsJac[iTerm].sum();
            else 
				ssdTermsJac[iTerm].x[0] += 1-ssdTermsJac[iTerm].sum();	

			// Stopping criteria
			double criteria = abs(QTerms.x[iTerm]-QTermPrev)/abs(QTerms.x[iTerm]);
			if(criteria<stopping){
				CONTINUE = false;
				cout<<iTerm<<"th term of Q has converged with "<<iIter<<" iterations\n";
			}
			QTermPrev = QTerms.x[iTerm];
				
		} // end of iterations
		cout<<"\nQTerm"<<"["<<iTerm<<"] = "<<QTerms.x[iTerm];
		gettimeofday(&dw2, 0);
		cout<<"\nTime for term is "<<(dw2.tv_sec - dw1.tv_sec)<<" seconds"<<endl;	
	} // finish calculating terms
	saveterms(fOutput, QTerms);
	cout<<"\nIt's over.";
//	getchar();
	return 0;
}
		
	
//			for(int iState = nStates-1; iState>0; iState--){
//				ssdTermsJac[iTerm].x[iState] = (1-W)*ssdTmp.x[iState] +W*ssdTermsJac[iTerm].x[iState];
//			}	
