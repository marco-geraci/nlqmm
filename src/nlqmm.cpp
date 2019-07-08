#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double C_hfun_nlqmm(NumericVector res, NumericVector u, NumericMatrix H, int N, int Q, double tau, double omicron){

	double fidelity = 0;
	double penalty = 0;
	double a = (tau - 1)*omicron;
	double b = tau*omicron;
	double val = 0;
	double Ajj = 0;
	double bvec = 0;
	double cvec = 0;
	
	// penalty
	
	for(int j = 0; j < Q; ++j){
		for(int k = 0; k < Q; ++k){
			penalty += u[k]*H(k,j)*u[j];
		}
	}

	// fidelity
	
	for(int i = 0; i < N; ++i){
		int s = 0;
		if(res[i] <= a){
			s = -1;
		} else if(res[i] >= b){
			s = 1;
		}
		Ajj = (1 - pow(s, 2))/omicron;
		bvec = s*((2*tau - 1)*s + 1);
		cvec = 0.5*(1-2*tau)*omicron*s - 0.5*(1-2*tau+2*pow(tau,2))*omicron*pow(s,2);
		
		fidelity += pow(res[i], 2)*Ajj + bvec*(res[i]) + cvec;
	}

	val = fidelity + penalty;
	return val;
}

// [[Rcpp::export]]
List C_loss_nlqmm(NumericVector res, NumericMatrix u, NumericVector weights, NumericMatrix Psiinv, double detH, double detPsi, int M, int N, int Q, double sigma, double tau, double omicron){

	NumericVector pen(M);
	double fidelity = 0;
	double penalty = 0;
	double a = (tau - 1)*omicron;
	double b = tau*omicron;
	double val = 0;
	double Ajj = 0;
	double bvec = 0;
	double cvec = 0;
	
	// penalty
	
	for(int i = 0; i < M; ++i){

		NumericVector tmp(Q);
		for(int j = 0; j < Q; ++j){
			for(int k = 0; k < Q; ++k){
				tmp[j] += u(i,k)*Psiinv(k,j)*u(i,j);
			}
			pen[i] += tmp[j];
		}
	penalty += pen[i];
	}
	
	// fidelity
	
	for(int i = 0; i < N; ++i){
		int s = 0;
		if(res[i] <= a){
			s = -1;
		} else if(res[i] >= b){
			s = 1;
		}
		Ajj = (1 - pow(s, 2))/omicron;
		bvec = s*((2*tau - 1)*s + 1);
		cvec = 0.5*(1-2*tau)*omicron*s - 0.5*(1-2*tau+2*pow(tau,2))*omicron*pow(s,2) - weights[i];
		
		fidelity += pow(res[i], 2)*Ajj + bvec*(res[i]) + cvec;
	}


	val = N * log(tau * (1- tau)/sigma) - 0.5*(M*detPsi + detH + (fidelity + penalty)/sigma);
	
	List ans;
	ans["val"] = val;
	ans["hsum"] = fidelity + penalty;
	return ans;
}

