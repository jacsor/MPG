#include "RcppArmadillo.h"
#include "helpers.h"
#include "gibbs.h"



using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
Rcpp::List fit( arma::mat Y, 
                arma::vec C, 
                Rcpp::List prior,
                Rcpp::List mcmc,
	              Rcpp::List state )
{
  Rcpp::RNGScope scope;  
  MCMC H(    Y, 
             C, 
             prior,
             mcmc,
             state  
	      );
        
  List chain = H.get_chain();      
        
  List data = Rcpp::List::create(  
      Rcpp::Named( "Y" ) = Y,
      Rcpp::Named( "C" ) = C
    ) ; 
  
  return Rcpp::List::create(  
      Rcpp::Named( "chain" ) = chain,
      Rcpp::Named( "data" ) = data,
      Rcpp::Named( "prior" ) = prior,
      Rcpp::Named( "mcmc" ) = mcmc
    ) ;    
}

