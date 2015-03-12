#ifndef HELPERS_H
#define HELPERS_H

#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;


double log_exp_x_plus_exp_y(double x, double y);

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma); 

arma::mat rWishartArma(arma::mat Sigma, int df);

arma::vec dmvnrm_arma_precision(  arma::mat x,  
                                  arma::rowvec mean,  
                                  arma::mat omega, 
                                  bool logd = true);
                                  
double rgammaBayes(double shape, double rate);


double beta_fun(arma::vec alpha, bool logB = true);

double marginalLikeDirichlet(arma::vec data, arma::vec alpha, bool logM = true);

arma::vec rDirichlet(arma::vec alpha, bool logR = true);

double dGeneralizedBeta(double x, double a, double b, arma::vec extremes, bool logD = true  );

int sampling(vec probs);


double KL(  arma::vec mu_1, 
            arma::vec mu_2, 
            arma::mat Sigma_1, 
            arma::mat Omega_1, 
            arma::mat Sigma_2, 
            arma::mat Omega_2);


#endif