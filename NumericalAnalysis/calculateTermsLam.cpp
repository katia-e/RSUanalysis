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
#include <random>


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
	double LamJac = 1.5;					// reference point for series expansion	
	double W = 0.9; 					// Stabilizing parameter
	
	for(int i = 1; i<argc; i++) 
		switch (i){			
			case 1: LamJac = atof(argv[i]);
			case 2: car_rate = atof(argv[i]);
			case 3: K = atoi(argv[i]); 
			case 4: C = atoi(argv[i]); 
			case 5: M = atoi(argv[i]); 
			case 6: nTerms = atoi(argv[i]);
			case 7: stopping = atof(argv[i]);
			case 8: W = atof(argv[i]);
			case 9: beta[0] = atof(argv[i]);
			case 10: beta[1] = atof(argv[i]);
			case 11: beta[2] = atof(argv[i]);
			case 12: data_rate[0] = atof(argv[i]);
			case 13: data_rate[1] = atof(argv[i]);
			case 14: data_rate[2] = atof(argv[i]);	
		};
	double percision = 1;
	string fOutput = "LamTermsK"+to_string(K)+
					"-M"+to_string(M)+"-C"+to_string(C)+
					"-lam"+toStringPrecision(LamJac,percision)+"-carRate" + toStringPrecision(car_rate,percision)+
					"-W"+toStringPrecision(W,percision) + "-stop"+toStringScient(stopping,0)+".dat";					
	cout<<endl<<"Lam0 = "<<LamJac<<endl
				<<"K = "<<K<<endl
				<<"C = "<<C<<endl
				<<"M = "<<M<<endl
				<<"Car Rate = "<<car_rate<<endl
				<<"Beta = "<<beta[0]<<", "<<beta[1]<<", "<<beta[2]<<endl
				<<"nTerms = "<<nTerms<<endl
				<<"eps = "<<stopping<<endl
				<<"W = "<<W<<endl
				<<"Output file is "<<fOutput<<endl;				


	/* Parameters of the state space */
	int DIM = M*(C+1);				// Number of dimentions in the queueing model
	cout<<"\nDimentions "<<DIM;	
	cout<<"\nDimentions "<<DIM;
	
	//--------------------------
	// Create state space
	int nStates = binomial(K, DIM);
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
	cout<<"\n Scheduler: ";
	print(mu);	
	
//	getchar();		
	vect *ssdTermsJac; ssdTermsJac = new vect [nTerms];	// Predefinitions for calculation
	vect QTerms; QTerms.init(nTerms);					// Terms of mean buffer content
	vect nCarsTerms; nCarsTerms.init(nTerms);			// Terms of mean number of cars in the system
	vect pAccTerms; pAccTerms.init(nTerms);				// Terms of accepted for transmission packets
	vect pTrTerms; pTrTerms.init(nTerms);				// Terms of rate for actually transmitted packets
	//======================================================================================================
	std::mt19937_64 rng;
	std::uniform_real_distribution<double> unif(0, 1.0/nStates);
	
	for(int iTerm=0; iTerm<nTerms; iTerm++){			// Initialisation of TERMS
		ssdTermsJac[iTerm].init(nStates);
		for(int iState=0; iState<nStates; iState++){
			ssdTermsJac[iTerm].x[iState] = 0;// nStates;//unif(rng);	// Initial guess
		}
	}
	ssdTermsJac[0].x[0] = 1;
	

	stState state, jstate; 
	state.init(DIM); jstate.init(DIM);
	timeval dw1, dw2;
	double creteria_previous = 0;
	for(int iTerm = 0; iTerm < nTerms; iTerm++){// start calculating terms
		gettimeofday(&dw1, 0);
		bool it_is_not_first_term = (iTerm>0);
		cout<<"\n"<<iTerm<<"th term";
		bool CONTINUE = true;
		int iIter = 0;
		double QTermPrev = 0;
		while(CONTINUE){
			iIter++;		
			bool THIS_IS_NOT_FIRST_ITER = true;
			if (iIter == 1)
				THIS_IS_NOT_FIRST_ITER = false;
			QTerms.x[iTerm] = 0;
			pAccTerms.x[iTerm] = 0;
			pTrTerms.x[iTerm] = 0;			
//		for(int iIter = 1; iIter<nIter; iIter++){// start Jacobi iterations
			//bool it_is_last_iteration = (iIter==nIter-1);
			state.to0(); state.st[0] = K;
			// ======= FORMULA IMPLEMENTATION pi[i](Q0+lam0*Q1) = -pi[i-1]Q1
			for(int iState = nStates-1; iState>0; iState--){ //calculation for each state
//			for(int iState = 0; iState< nStates-1; iState++){ //calculation for each state
				int nCars = state.sum();
				double Pr = 0; double bi = 0; double aii = 0;
				double tmp = 0; int x; int from, to,up, down; // sup variable
				int jState;
				//==========================================================
				// Solution for Ax = b
				//  All possible transissions from iState a11 - Diagonal of A = Q0+lam0*Q1
				if(nCars < K){ // Get from iState by a car arrival
					aii -= car_rate;
				}
				// Get from iState by a packet arrival or car range switch
				for(int n = 0; n < C+1; n++)					
					for(int m = 0; m < M; m++){
						x = state.st[M*n+m];
						if(n<C)  
							aii -= x*(LamJac+beta[m]);
						else           
							aii -= x*beta[m];
						}
				// Get from iState by a packet departure
				if(state.sum(M,DIM)>0){
					tmp = 0; 
					for(int n=1;n<C+1;n++)
						for(int m=0;m<M;m++){
							int ind = M*n+m;
							tmp += state.st[ind]*mu.x[ind-M];
						}
					aii -= tmp;
				}
				//============================================================
				// Get to iState by a car arrival
				if(state.st[0]>0){ // then getting to iState by a car arrival is possible
					jstate.copy(state);
					jstate.st[0] --; // % From this state	
					jState = getStateNum(StateSpace, jstate);				
					Pr += car_rate*ssdTermsJac[iTerm].x[jState];
				}			
				if (it_is_not_first_term){// Perturbed component 1
					for(int n = 0; n < C; n++)					
						for(int m = 0; m < M; m++){
							x = state.st[M*n+m];  
							bi += x*ssdTermsJac[iTerm-1].x[iState];	
						}			
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
								tmp+=LamJac* x* ssdTermsJac[iTerm].x[jState];
                                if(it_is_not_first_term) // Perturbed part component 2
                                    bi -= x*ssdTermsJac[iTerm-1].x[jState]; 								
							}
						}
					}
					Pr += tmp;
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
								tmp += x*mu_scheduled*ssdTermsJac[iTerm].x[jState];//flag
							}
						}
					 }
					Pr = Pr + tmp; 
				}			       
//				ssdTmp.x[iState] = (bi-Pr)/aii;; // Store previous value
//				ssdTermsJac[iTerm].x[iState] = (bi-Pr)/aii;  
				// =============== Weighted Jacobi method ========================
				ssdTermsJac[iTerm].x[iState] =  (1-W)*ssdTermsJac[iTerm].x[iState]+W*(bi-Pr)/aii;	
				//================================================================
				//================= PERFORMANCE MEASURES =========================
				int nPacketsInStates = 0;				
				int rateTrPackets = 0;
				for(int n = 1; n<C+1; n++)
					nPacketsInStates += state.sum(n*M,(n+1)*M)*n;						
				for(int m = 0; m<M; m++)
					for(int n = 1; n<C+1; n++){
						int position = M*n+m;
						rateTrPackets += state.st[position]*mu.x[position-M];
					}
				QTerms.x[iTerm] += nPacketsInStates*ssdTermsJac[iTerm].x[iState];					
				nCarsTerms.x[iTerm] += state.sum()*Pr;
				pAccTerms.x[iTerm] 	+= ssdTermsJac[iTerm].x[iState]*state.sum(0,DIM-M);
				pTrTerms.x[iTerm] 	+= ssdTermsJac[iTerm].x[iState]*rateTrPackets;
//				incrvect(state, state, K);
				dcrvect(state, state, K);
			} // all states are calculated
			// normalize terms
            if(it_is_not_first_term) // sum(term) = 0
                ssdTermsJac[iTerm].x[0] -= ssdTermsJac[iTerm].sum();
            else 
				ssdTermsJac[iTerm].x[0] += 1-ssdTermsJac[iTerm].sum();	

			// Stopping criteria			
			double  criteria = abs(QTerms.x[iTerm]-QTermPrev)/abs(QTerms.x[iTerm]);
//			if ((criteria>creteria_previous) && THIS_IS_NOT_FIRST_ITER){
//				cout<<"doesn't converge at all "<<criteria<<endl;
//				getchar();
//			}
			creteria_previous = criteria;
			if (QTerms.x[iTerm]==0){
				cout<<endl;
				CONTINUE = false;
			}
			
//			cout<<criteria<<endl;
			if(criteria<stopping){
				CONTINUE = false;
				cout<<iTerm<<"th term of Q has converged with "<<iIter<<" iterations\n";
				cout<<"QTerm["<<iTerm<<"] = "<<QTerms.x[iTerm]<<endl;
				cout<<iTerm<<"th term of rate of accepted packets"<<endl;
				cout<<"pAccTerms["<<iTerm<<"] = "<<pAccTerms.x[iTerm]<<endl;
				cout<<iTerm<<"th term of rate of transmitted packets"<<endl;
				cout<<"pAccTerms["<<iTerm<<"] = "<<pTrTerms.x[iTerm]<<endl;				
				}
			QTermPrev = QTerms.x[iTerm];

				
		} // end of iterations
//		cout<<"\nQTerm"<<"["<<iTerm<<"] = "<<QTerms.x[iTerm];
		gettimeofday(&dw2, 0);
		cout<<"\nTime for term is "<<(dw2.tv_sec - dw1.tv_sec)<<" seconds"<<endl;	
	} // finish calculating terms
	saveterms(fOutput, QTerms,pAccTerms,pTrTerms);
	cout<<"\nIt's over.";
//	getchar();
	return 0;
}
		
	
//			for(int iState = nStates-1; iState>0; iState--){
//				ssdTermsJac[iTerm].x[iState] = (1-W)*ssdTmp.x[iState] +W*ssdTermsJac[iTerm].x[iState];
//			}	
