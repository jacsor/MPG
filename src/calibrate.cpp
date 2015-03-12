#include "helpers.h"
#include "RcppArmadillo.h"


using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List calib( arma::mat Y, 
                  arma::vec C,
                  arma::mat Z,  
                  NumericVector mu_input,
                  IntegerVector mu_dim,
                  NumericVector mu0_input,
                  IntegerVector mu0_dim,
                  arma::mat S )
{
  
  arma::cube mu(mu_input.begin(), mu_dim[0], mu_dim[1], mu_dim[2]);
  arma::cube mu0(mu0_input.begin(), mu0_dim[0], mu0_dim[1], mu0_dim[2]);
  
  int n = Y.n_rows;
  int p = Y.n_cols;
  int niter = Z.n_rows;
  cube calibration(niter,p,n); calibration.fill(0); 
  mat calibrationMedian(n,p); calibrationMedian.fill(0);
  vec temp(n);
  
  for(int it=0; it<niter; it++)
  {
    for(int i=0; i<n; i++)
    {
      if(S(it,Z(it,i)) == 1)
      {
        calibration.slice(i).row(it) =   mu.slice(it).cols(Z(it,i)*p,Z(it,i)*p+p-1).row(C(i)) - mu0.slice(it).col(Z(it,i)).t();
        //calibrationMean.row(i) += calibration.slice(it).row(i)/niter;
      }
    }    
  }  
  
  for( int i = 0; i < n; i++)
  {
    calibrationMedian.row(i) = median(calibration.slice(i),0);
  }
  
  
  return Rcpp::List::create(  
    Rcpp::Named( "Y_cal" ) = Y - calibrationMedian,
    Rcpp::Named( "calibration_distribution" ) = calibration,
    Rcpp::Named( "calibration_median" ) = calibrationMedian 
  ) ;    
}
