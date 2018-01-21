#include <fstream>
#include <iomanip> // set precision
#include "math.h"
#include "stdint.h"
const int PRECISION = 5;

struct stState{
	uint16_t *st;
	uint16_t n;
	void init(int N){
		n = N;
		st = new uint16_t [n];
		for(int i = 0; i<n; i++)
			st[i] = 0;
	}
	void copy(stState a){
		// check
		if (n < a.n){
			printf ("stState: Vector should have the same or smaller length than the storage");
			exit (EXIT_FAILURE);
		  }
		for(int i = 0; i<n; i++){
			n = a.n;
			st[i] = a.st[i];
		}
	}
	void to0(){
		for(int i = 0; i<n; i++)
			st[i] = 0;
	}
	int sum(){
		int res = 0;
		for(int i = 0; i<n; i++)
			res += st[i];
		return res;
	}
	int sum(int a, int b){
		if (n < b){
			printf ("stState: intex is out of range");
			exit (EXIT_FAILURE);
		  }
		int res = 0;
		for(int i = a; i<b; i++)
			res += st[i];
		return res;
	}
	stState substate(int N){
		stState res;
		res.n = N;
		res.st = st;
		return res;
	};
	uint16_t numOfnon0(){
		uint16_t res = 0;
		for(int i = 0; i<n; i++ )
			if(st[i]>0)
				res++;
		return res;
	}
};

struct vect{
	double *x;
	int n;
	void init(int N){
		n = N;
		x = new double [n];
		for(int i = 0; i<n; i++)
			x[i] = 0;
	}
	
	void fill(double *sub_vector, int size){
	// Fill values of @x with values of given array @sub_vector of size @size
	// $n must be devidible by @size
	if (n % size!= 0) {
		cout<<"ERR vect.fill: $n must be devidible by @size";
		exit (EXIT_FAILURE);
	}
	for (int i=0; i<n/size; i++)
		for (int j=0; j<size; j++)
			x[i*size+j] = sub_vector[j];
	}
	void copy(vect a){
		if (n < a.n){
			printf ("stState: Vector being copied should have the same or smaller length than the storage");
			exit (EXIT_FAILURE);
		  }
		for(int i = 0; i<n; i++){
			n = a.n;
			x[i] = a.x[i];
		}
	}
	void to0(){
		for(int i = 0; i<n; i++)
			x[i] = 0;
	}
	double sum(){
		double res = 0;
		for(int i = 0; i<n; i++)
			res += x[i];
		return res;
	}
	vect substate(int N){
		vect res;
		res.n = N;
		res.x = x;
		return res;
	};
	void print(){
		cout<<"\n";
		for(int i = 0; i<n; i++)
			cout<<x[i]<<"\t";		
	}
	void norm(){// forse sum to 1
		double s; s = sum();
		for(int i = 0; i<n; i++)
			x[i] /= s;		
	}
	void mult(double a){
		for(int i = 0; i<n; i++)
			x[i] *= a;
	}
};
int binomial(int K, int DIM){ // Calculates the number of states in the state space
	int res = K+1;
	for(int i = 2; i <= DIM; i++)
		res = res*(K+i)/i; 
	return res;
}

int convertTo10(int num, int base){ // convert between 10 and base numerical systems
	int num10 = 0; bool DO = true; int res = 0; int i = 0;
	while (num){		 
		if (num>=10) res = num%10;
		else res = num;
		num10 += res*pow(base,i);
		num = num/10; i++;
	}
	return num10;
}

void incrvect(stState nextState, stState state, int K){ // Finding the next state of the system
//	int base = K+1; int res = 0;
	nextState.copy(state);
	int nLast = state.n-1;
	if (nextState.sum() == K){
		incrvect(nextState, nextState.substate(nextState.n-1),K);
		nextState.st[nLast] = 0;
	}
	else{
		if (nextState.st[nLast]<K)
			nextState.st[nLast]++;
		else{
			incrvect(nextState, nextState.substate(nextState.n-1),K-1);
			nextState.st[nLast] = 0;
		}
	}
}

void dcrvect(stState nextState, stState state, int K){ // Finding the previous state of the system
	nextState.copy(state);
	int nLast = state.n-1;
	if (nextState.st[nLast] == 0){
		nextState.st[nLast] = K;
		dcrvect(nextState, nextState.substate(nextState.n-1), K);
		while(nextState.sum()>K)
			dcrvect(nextState, nextState, K);
	}
	else
		nextState.st[nLast]-=1;
}

void linspace(vect mu, double a, double b){
	int n = mu.n;
	if(n==1)
		mu.x[0] = a;
	else{
		double delta = (b-a)/(n-1);
		for(int i=0; i<n; i++)
			mu.x[i] = a+i*delta;
	}		
}
void scheduler00(vect mu, double muSim){
	linspace(mu, 0.5, 1);
	mu.norm();
	mu.mult(muSim);
}

void scheduler03(vect mu, int C){
	int M = 3;
	mu.x[0] = 0.1;
	mu.x[1] = 1;
	mu.x[2] = 0.1;
	for(int n = 1; n<C; n++)
		for(int m = 0; m<M; m++)
			mu.x[M*n+m] = mu.x[m]; 
}

void scheduler05(vect mu, int C){
	int M = 5;
	mu.x[0] = 0.5;
	mu.x[1] = 0.75;
	mu.x[2] = 1;
	mu.x[3] = 0.75;
	mu.x[4] = 0.5;	
	for(int n = 1; n<C; n++)
		for(int m = 0; m<M; m++)
			mu.x[M*n+m] = mu.x[m]; 	
}

void saveterms(string f_name, vect terms){	
	ofstream file;
	file.open(f_name);
	int nTerms = terms.n;
	for(int iTerm = 0; iTerm < nTerms; iTerm++)
		file<<terms.x[iTerm]<<std::setprecision(PRECISION)<<"\n";	
	file.close();
}

// 12 Dimentions
void allocateStateSpace(int************X, int K){
	int DIM = 12; int iState = 0;
	stState ind; ind.init(DIM);
	for(ind.st[0] = 0; ind.st[0] < K+1; ind.st[0]++){
		int size = K-ind.st[0]+1;
		X[ind.st[0]] = new int********** [size];
		for(ind.st[1] = 0; ind.sum(0,2) < K+1; ind.st[1]++){
			size = K-ind.sum(0,2)+1;
			X[ind.st[0]][ind.st[1]] = new int********* [size];
			for(ind.st[2] = 0; ind.sum(0,3) < K+1; ind.st[2]++){
				size = K-ind.sum(0,3)+1;
				X[ind.st[0]][ind.st[1]][ind.st[2]] = new int******** [size];
				for(ind.st[3] = 0; ind.sum(0,4) < K+1; ind.st[3]++){
					size = K-ind.sum(0,4)+1;
					X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]] = new int******* [size];
					for(ind.st[4] = 0; ind.sum(0,5) < K+1; ind.st[4]++){
						size = K-ind.sum(0,5)+1;
						X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]] = new int****** [size];
						for(ind.st[5] = 0; ind.sum(0,6) < K+1; ind.st[5]++){
							size = K-ind.sum(0,6)+1;
							X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]] = new int***** [size];
							for(ind.st[6] = 0; ind.sum(0,7) < K+1; ind.st[6]++){
								size = K-ind.sum(0,7)+1;
								X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]][ind.st[6]] = new int**** [size];
								for(ind.st[7] = 0; ind.sum(0,8) < K+1; ind.st[7]++){
									size = K-ind.sum(0,8)+1;
									X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]][ind.st[6]][ind.st[7]] = new int*** [size];
									for(ind.st[8] = 0; ind.sum(0,9) < K+1; ind.st[8]++){
										size = K-ind.sum(0,9)+1;
										X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]][ind.st[6]][ind.st[7]][ind.st[8]] = new int** [size];
										for(ind.st[9] = 0; ind.sum(0,10) < K+1; ind.st[9]++){
											size = K-ind.sum(0,10)+1;
											X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]][ind.st[6]][ind.st[7]][ind.st[8]][ind.st[9]] = new int* [size];
											for(ind.st[10] = 0; ind.sum(0,11) < K+1; ind.st[10]++){
												size = K-ind.sum(0,11)+1;
												X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]][ind.st[6]][ind.st[7]][ind.st[8]][ind.st[9]][ind.st[10]] = new int [size];
												for(ind.st[11] = 0; ind.sum()< K+1; ind.st[11]++){
													X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]][ind.st[6]][ind.st[7]][ind.st[8]][ind.st[9]][ind.st[10]][ind.st[11]] = iState;
													iState++;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

// 10 Dimentions
void allocateStateSpace(int**********X, int K){
	int DIM = 12; int iState = 0;
	stState ind; ind.init(DIM);
	for(ind.st[0] = 0; ind.st[0] < K+1; ind.st[0]++){
		int size = K-ind.st[0]+1;
		X[ind.st[0]] = new int******** [size];
		for(ind.st[1] = 0; ind.sum(0,2) < K+1; ind.st[1]++){
			size = K-ind.sum(0,2)+1;
			X[ind.st[0]][ind.st[1]] = new int******* [size];
			for(ind.st[2] = 0; ind.sum(0,3) < K+1; ind.st[2]++){
				size = K-ind.sum(0,3)+1;
				X[ind.st[0]][ind.st[1]][ind.st[2]] = new int****** [size];
				for(ind.st[3] = 0; ind.sum(0,4) < K+1; ind.st[3]++){
					size = K-ind.sum(0,4)+1;
					X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]] = new int***** [size];
					for(ind.st[4] = 0; ind.sum(0,5) < K+1; ind.st[4]++){
						size = K-ind.sum(0,5)+1;
						X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]] = new int**** [size];
						for(ind.st[5] = 0; ind.sum(0,6) < K+1; ind.st[5]++){
							size = K-ind.sum(0,6)+1;
							X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]] = new int*** [size];
							for(ind.st[6] = 0; ind.sum(0,7) < K+1; ind.st[6]++){
								size = K-ind.sum(0,7)+1;
								X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]][ind.st[6]] = new int** [size];
								for(ind.st[7] = 0; ind.sum(0,8) < K+1; ind.st[7]++){
									size = K-ind.sum(0,8)+1;
									X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]][ind.st[6]][ind.st[7]] = new int* [size];
									for(ind.st[8] = 0; ind.sum(0,9) < K+1; ind.st[8]++){
										size = K-ind.sum(0,9)+1;
										X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]][ind.st[6]][ind.st[7]][ind.st[8]] = new int [size];
										for(ind.st[9] = 0; ind.sum(0,10) < K+1; ind.st[9]++){
											X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]][ind.st[6]][ind.st[7]][ind.st[8]][ind.st[9]]= iState;
											iState++;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
void allocateStateSpace(int******X, int K){ // for example C = 2; M = 2
	int DIM = 6; int iState = 0;
	stState ind; ind.init(DIM);
	for(ind.st[0] = 0; ind.st[0] < K+1; ind.st[0]++){
		int size = K-ind.st[0]+1;
		X[ind.st[0]] = new int**** [size];
		for(ind.st[1] = 0; ind.sum(0,2) < K+1; ind.st[1]++){
			size = K-ind.sum(0,2)+1;
			X[ind.st[0]][ind.st[1]] = new int*** [size];
			for(ind.st[2] = 0; ind.sum(0,3) < K+1; ind.st[2]++){
				size = K-ind.sum(0,3)+1;
				X[ind.st[0]][ind.st[1]][ind.st[2]] = new int** [size];
				for(ind.st[3] = 0; ind.sum(0,4) < K+1; ind.st[3]++){
					size = K-ind.sum(0,4)+1;
					X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]] = new int* [size];
					for(ind.st[4] = 0; ind.sum(0,5) < K+1; ind.st[4]++){
						size = K-ind.sum(0,5)+1;
						X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]] = new int [size];
						for(ind.st[5] = 0; ind.sum(0,6) < K+1; ind.st[5]++){
							X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]] = iState;
							iState++;
						}
					}
				}
			}
		}
	}
}

int getStateNum(int************X, stState ind){ //12 Dimentions
	return X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]][ind.st[6]][ind.st[7]][ind.st[8]][ind.st[9]][ind.st[10]][ind.st[11]];
}

int getStateNum(int**********X, stState ind){// 10 Dimentions
	return X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]][ind.st[6]][ind.st[7]][ind.st[8]][ind.st[9]];
}

int getStateNum(int******X, stState ind){// 6 Dimentions
	return X[ind.st[0]][ind.st[1]][ind.st[2]][ind.st[3]][ind.st[4]][ind.st[5]];
}

void print(vect x){
	int len = x.n;
	for(int i=0; i<len; i++)
		cout<<x.x[i]<<"\t";
}

void print(vect* ssdTerms, int nTerms){
	int len = ssdTerms[0].n;
	for(int iState=0; iState<len; iState++){
		for(int iTerm=0; iTerm<nTerms; iTerm++)
			cout<<ssdTerms[iTerm].x[iState]<<"\t";
		cout<<"\n";
	}
}
