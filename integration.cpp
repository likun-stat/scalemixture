//
//  main.cpp
//  integration
//
//  Created by ZhangLikun on 9/17/18.
//  Copyright © 2018 ZhangLikun. All rights reserved.
//



#include <RcppGSL.h>
#include <math.h>
#include <gsl/gsl_integration.h>

// [[Rcpp::depends(RcppGSL)]]

using namespace Rcpp;


struct my_f_params {double x; double tau_sqd; double delta;};


// xval > 1
double mix_me_distn_integrand(double x, void *p) {
    struct my_f_params *params = (struct my_f_params *)p;
    double xval   = (params->x);
    double tau_sqd  = (params->tau_sqd);
    double delta = (params->delta);
    
    double half_result = (delta/(2*delta-1))*pow(xval-x,-(1-delta)/delta)-((1-delta)/(2*delta-1))*pow(xval-x,-1);
    
    return R::dnorm(x, 0.0, sqrt(tau_sqd), 0) * half_result;
}


double mix_me_dens_integrand(double x, void *p) {
    struct my_f_params *params = (struct my_f_params *)p;
    double xval   = (params->x);
    double tau_sqd  = (params->tau_sqd);
    double delta = (params->delta);
    
    double half_result = ((1-delta)/(2*delta-1))*pow(xval-x,-2)-((1-delta)/(2*delta-1))*pow(xval-x,-1/delta);

    return R::dnorm(x, 0.0, sqrt(tau_sqd), 0) * half_result;
}


// [[Rcpp::export]]
double pmixture_me_uni(double x, double tau_sqd, double delta, double relerr = 1e-10) {
    
    double result = 0.0;
    
    gsl_function F;
    F.function = &mix_me_distn_integrand;
    
    gsl_integration_workspace *work = gsl_integration_workspace_alloc(1e5);
    
    gsl_set_error_handler_off();
    
    struct my_f_params params = { x, tau_sqd, delta };
    F.params = &params;
    
    double abserr = 0.0;
    
    if(x>1){
        
        // QAGI adaptive integration on infinite intervals
        double err = gsl_integration_qagil(&F, x-1, 1e-12, relerr, 1e5, work, &result, &abserr);
        
        if (!ISNAN(err)){
            result = R::pnorm(x-1, 0.0, sqrt(tau_sqd), 1, 0) - result;
        }
        else {
            Rcpp::Rcout << "Error in integration. Returning -1" << std::endl;
            Rcpp::Rcout << "Err = " << err << std::endl;
            result = -1.0;
        }
    }
    else {result = -1.0;}
    
    gsl_integration_workspace_free(work);
    
    return result;
}



// [[Rcpp::export]]
NumericVector pmixture_me_old(NumericVector x, double tau_sqd, double delta, double relerr = 1e-10) {
    
    int n = x.size();
    NumericVector resultVec(n);
    
    gsl_function F;
    F.function = &mix_me_distn_integrand;
    
    gsl_integration_workspace *work = gsl_integration_workspace_alloc(1e5);
    
    gsl_set_error_handler_off();
    
    for(int i = 0; i < n; i++) {
        struct my_f_params params = { x[i], tau_sqd, delta };
        F.params = &params;
        
        double result = 0.0;
        double abserr = 0.0;
        
        if(x[i]>1){
            
            // QAGI adaptive integration on infinite intervals
            double err = gsl_integration_qagil(&F, x[i]-1, 1e-12, relerr, 1e5, work, &result, &abserr);
            
            if (!ISNAN(err)){
                result = R::pnorm(x[i]-1, 0.0, sqrt(tau_sqd), 1, 0) - result;
            }
            else {
                Rcpp::Rcout << "Error in integration. Returning -1" << std::endl;
                Rcpp::Rcout << "Err = " << err << std::endl;
                result = -1.0;
            }
        }
        else {result = -1.0;}
        
        resultVec[i] = result;
    }
    
    gsl_integration_workspace_free(work);
    
    return resultVec;
}


// [[Rcpp::export]]
double asymptotic_p(double x, double delta){
    double result = 1-(delta/(2*delta-1))*pow(x,(delta-1)/delta)+((1-delta)/(2*delta-1))*pow(x,-1);
    return result;
}


// [[Rcpp::export]]
NumericVector pmixture_me(NumericVector x, double tau_sqd, double delta, double relerr = 1e-10) {
    
    int n = x.size();
    NumericVector resultVec(n);
    
    gsl_function F;
    F.function = &mix_me_distn_integrand;
    
    gsl_integration_workspace *work = gsl_integration_workspace_alloc(1e5);
    
    gsl_set_error_handler_off();
    
    
    // First decide a cutoff point where the quantile is high enough to use asymptotic results
    double high_quantiles[] = {10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 15.0, 20.0, 50.0};
    double threshold = 50.0;
    double init_p = 0.0;
    double init_est = 0.0;
    
    for(int iter = 0; iter < 12; iter++){
        init_p = pmixture_me_uni(high_quantiles[iter], tau_sqd, delta);
        init_est = asymptotic_p(high_quantiles[iter], delta);
        if(fabs(init_p - init_est)<0.0005){
            threshold = high_quantiles[iter];
            break;
        }
    }
    
    // Calculate CDF fucntion at each x values
    for(int i = 0; i < n; i++) {
        struct my_f_params params = { x[i], tau_sqd, delta };
        F.params = &params;
        
        double result = 0.0;
        double abserr = 0.0;
        
        if(x[i]>1 & x[i]<threshold){
            
            // QAGI adaptive integration on infinite intervals
            double err = gsl_integration_qagil(&F, x[i]-1, 1e-12, relerr, 1e5, work, &result, &abserr);
            
            if (!ISNAN(err)){
                result = R::pnorm(x[i]-1, 0.0, sqrt(tau_sqd), 1, 0) - result;
            }
            else {
                Rcpp::Rcout << "Error in integration. Returning -1" << std::endl;
                Rcpp::Rcout << "Err = " << err << std::endl;
                result = -1.0;
            }
        }
        else if(x[i]>=threshold) {result = asymptotic_p(x[i], delta);}
        else {result = -1.0;}
        
        resultVec[i] = result;
    }
    
    gsl_integration_workspace_free(work);
    
    return resultVec;
}


// [[Rcpp::export]]
NumericVector dmixture_me_old(NumericVector x, double tau_sqd, double delta, double relerr = 1e-6) {
    
    int n = x.size();
    NumericVector resultVec(n);
    
    gsl_function F;
    F.function = &mix_me_dens_integrand;
    
    gsl_integration_workspace *work = gsl_integration_workspace_alloc(1e5);
    
    gsl_set_error_handler_off();
    
    for(int i = 0; i < n; i++) {
        struct my_f_params params = { x[i], tau_sqd, delta };
        F.params = &params;
        
        double result = 0.0;
        double abserr = 0.0;
        
        if(x[i]>1){
            
            // QAGI adaptive integration on infinite intervals
            double err = gsl_integration_qagil(&F, x[i]-1, 1e-12, relerr, 1e5, work, &result, &abserr);
            
            if (!ISNAN(err)){
                result = - result;
            }
            else {
                Rcpp::Rcout << "Error in integration. Returning -1" << std::endl;
                Rcpp::Rcout << "Err = " << err << std::endl;
                result = -1.0;
            }
        }
        else {result = -1.0;}
        
        if (result > 0) {
            resultVec[i] = result;
        } else {
            resultVec[i] = 0;
        }
    }
    
    gsl_integration_workspace_free(work);
    
    return resultVec;
}



// [[Rcpp::export]]
double asymptotic_d(double x, double delta){
    double result = ((1-delta)/(2*delta-1))*(pow(x,-1/delta)-pow(x,-2));
    return result;
}



// [[Rcpp::export]]
NumericVector dmixture_me(NumericVector x, double tau_sqd, double delta, double relerr = 1e-6) {
    
    int n = x.size();
    NumericVector resultVec(n);
    
    gsl_function F;
    F.function = &mix_me_dens_integrand;
    
    gsl_integration_workspace *work = gsl_integration_workspace_alloc(1e5);
    
    gsl_set_error_handler_off();
    
    // First decide a cutoff point where the quantile is high enough to use asymptotic results
    double high_quantiles[] = {10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 15.0, 20.0, 50.0};
    double threshold = 50.0;
    double init_p = 0.0;
    double init_est = 0.0;
    
    for(int iter = 0; iter < 12; iter++){
        init_p = pmixture_me_uni(high_quantiles[iter], tau_sqd, delta);
        init_est = asymptotic_p(high_quantiles[iter], delta);
        if(fabs(init_p - init_est)<0.0005){
            threshold = high_quantiles[iter];
            break;
        }
    }
    
    // Calculate PDF fucntion at each x values
    for(int i = 0; i < n; i++) {
        struct my_f_params params = { x[i], tau_sqd, delta };
        F.params = &params;
        
        double result = 0.0;
        double abserr = 0.0;
        
        if(x[i]>1 & x[i]<threshold){
            
            // QAGI adaptive integration on infinite intervals
            double err = gsl_integration_qagil(&F, x[i]-1, 1e-12, relerr, 1e5, work, &result, &abserr);
            
            if (!ISNAN(err)){
                result = - result;
            }
            else {
                Rcpp::Rcout << "Error in integration. Returning -1" << std::endl;
                Rcpp::Rcout << "Err = " << err << std::endl;
                result = -1.0;
            }
        }
        else if (x[i]>threshold) {result = asymptotic_d(x[i], delta);}
        else {result = -1.0;}
        
        if (result > 0) {
            resultVec[i] = result;
        } else {
            resultVec[i] = 0;
        }
    }

    gsl_integration_workspace_free(work);
    
    return resultVec;
}



/*  Test in R
 # parameter settings
 delta <- 0.3
 tau <- 1   # Always variance = std^2
 
 x_vals <- seq(1.001, 20, length = 10000)
 system.time(pdf_vals<-dmixture_me(x_vals, tau, delta))
 system.time(cdf_vals<-pmixture_me(x_vals, tau, delta))
 
 par(mfrow=c(2,1))
 plot(x_vals, pdf_vals, type="l")
 plot(x_vals, cdf_vals, type="l")
 */






// -------------------------------------------------------------------------- //
// This only makes sense if we want to search over x > 0.  Since we are interested
// in extremes, this seems okay -- things are not so extreme if they are in the
// lower half of the support of the (copula) distribution.
//                                                                            //

// [[Rcpp::export]]
NumericVector find_xrange_pmixture_me(double min_p, double max_p,
                                      NumericVector x_init,
                                      double tau_sqd, double delta, double relerr = 1e-10) {
    NumericVector min_x(1);
    NumericVector max_x(1);
    double p_min_x;
    double p_max_x;
    NumericVector x_range(2);
    
    min_x[0] = x_init[0];
    max_x[0] = x_init[1];
    
    if((min_x[0] <= 0) || (min_p <= 0.5)) {
        Rcpp::stop("This will only work for x > 0, which corresponds to p > 1/2.");
    }
    if(min_x[0] >= max_x[0]) {
        Rcpp::stop("x_init[1] must be smaller than x_init[2].");
    }
    
    
    // First the min
    p_min_x = pmixture_me(min_x, tau_sqd, delta, relerr)[0];
    while (p_min_x > min_p) {
        // Rcpp::Rcout << "left x is " << min_x[0] << std::endl;
        // Rcpp::Rcout << "F(" << min_x[0] << ") = " << p_min_x << std::endl;
        min_x[0] = min_x[0]/2;
        p_min_x = pmixture_me(min_x, tau_sqd, delta, relerr)[0];
    }
    x_range[0] = min_x[0];
    
    // Now the max
    p_max_x = pmixture_me(max_x, tau_sqd, delta, relerr)[0];
    while (p_max_x < max_p) {
        // Rcpp::Rcout << "right x is " << max_x << std::endl;
        // Rcpp::Rcout << "F(" << max_x << ") = " << p_max_x<< std::endl;
        max_x[0] = max_x[0]*2;
        p_max_x = pmixture_me(max_x, tau_sqd, delta, relerr)[0];
    }
    x_range[1] = max_x[0];
    
    return x_range;
}
//                                                                            //
// -------------------------------------------------------------------------- //




