// -------------------------------------------------------------------------- //
// -------------------------     Matrix Algebra     ------------------------- //
// -------------------------------------------------------------------------- //
//  Created by ZhangLikun on 9/17/18.
//  Copyright © 2018 ZhangLikun. All rights reserved.
//




#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export()]]
double eig2logdet_c (arma::vec d){
    arma::vec log_d = log(d);
    return sum(log_d);
}

// [[Rcpp::export()]]
double eig2inv_quadform_vector (arma::mat V, arma::vec d_inv, arma::vec x) {
    arma::mat step_r = V.t()*x;
    arma::mat step_m = diagmat(d_inv)*step_r;
    return dot(step_r.t(), step_m);
}

// [[Rcpp::export()]]
arma::mat eig2inv_times_vector (arma::mat V, arma::vec d_inv, arma::vec x) {
    arma::mat step_r = V.t()*x;
    arma::mat step_m = diagmat(d_inv)*step_r;

    return V*step_m;
}

// [[Rcpp::export()]]
double dmvn_eig (arma::mat R, arma::mat V, arma::vec d_inv){
    int n_rep = R.n_cols;
    double left = -0.5*n_rep*eig2logdet_c(1/d_inv);
    
    arma::vec right(n_rep);
    for(int i=0; i<n_rep; i++){
        arma::vec R_i = R.col(i);
        right(i) = eig2inv_quadform_vector(V, d_inv, R_i);
    }
    
    return left-0.5* sum(right);
}

// [[Rcpp::export()]]
arma::vec trans_NVec (NumericVector x){
    arma::vec y = as<arma::vec>(x);
    return y;
}

// [[Rcpp::export()]]
double X_s_likelihood_conditional(arma::vec X_s, double R, arma::mat V, arma::vec d){
    NumericVector tmp = qnorm(as<NumericVector>(wrap(1-R/X_s)));
    arma::vec X_s_to_Z = trans_NVec(tmp);
    double part1 = -0.5*eig2inv_quadform_vector(V, 1/d, X_s_to_Z);
    double part2 = 0.5*sum(X_s_to_Z % X_s_to_Z)-2*sum(log(X_s));
    return part1+part2;
}



// [[Rcpp::export()]]
double X_s_likelihood_conditional_on_X (arma::vec X_s, arma::vec X, double R, arma::mat V, arma::vec d, double tau_sqd){
    
    // treat X_s and X as vec
    // arma::vec X_s = vectorise(X_s_m);
    // arma::vec X = vectorise(X_m);
    
    if(any(X_s<R)) return -std::numeric_limits<double>::infinity();
    else{
        NumericVector tmp = qnorm(as<NumericVector>(wrap(1-R/X_s)));
        arma::vec X_s_to_Z = trans_NVec(tmp);
        double part1 = -0.5*eig2inv_quadform_vector(V, 1/d, X_s_to_Z);
        double part2 = 0.5*sum(X_s_to_Z % X_s_to_Z)-2*sum(log(X_s));
        double part3 = - 0.5*sum((X_s-X) % (X_s-X))/tau_sqd;
        
        return part1+part2+part3;
    }
}



// [[Rcpp::export()]]
List var_at_a_time_update_X_s (arma::vec X, double R, arma::mat V, arma::vec d, double tau_sqd, double v_q=0.5, int n_chain=100){
    
    int n_s = X.n_elem;
    arma::vec X_s = X;
    arma::uvec tmp =find(X_s<R);
    
    if(tmp.n_elem>0){
        for(int i=0; i<tmp.n_elem; i++){
            int ind = tmp(i);
            arma::vec nugget = 0.5*arma::randn(1);
            X_s(ind) = R + fabs(as_scalar(nugget));
        }
    }
    
    arma::vec accept(n_s);
    accept.fill(0);
    
    for(int i=0; i<n_chain; i++){
        for(int iter=0; iter<n_s; iter++){
            arma::vec X_s_update = X_s;
            X_s_update(iter) = X_s[iter]+v_q*as_scalar(arma::randn(1));
            double log_num = X_s_likelihood_conditional_on_X(X_s_update, X, R, V, d, tau_sqd);
            double log_denom = X_s_likelihood_conditional_on_X(X_s, X, R, V, d, tau_sqd);
            
            double r = exp(log_num - log_denom);
            if(as_scalar(arma::randu(1))<r){
                X_s = X_s_update;
                accept(iter) = accept(iter) + 1;
            }
        }
    }
    
    List result;
    result["X.s"] = X_s;
    result["accept"] = accept;
    
    return result;
}



