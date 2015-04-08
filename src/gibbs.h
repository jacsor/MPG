#ifndef GIBBS_H
#define GIBBS_H

#include "RcppArmadillo.h"
#include <omp.h>

using namespace Rcpp;
using namespace arma;


class MCMC
{
  private:
  mat Y;            // the data (n,J)
  vec C;            // the group membership (n,1)
  int K_0, K_1;     // number of mixture components 
  int J;            // the number of groups
  int n;            // number of observations
  int p;            // observation dimension
  int num_iter, num_burnin, num_thin, num_display;   // number of iteration, burnin and thinning
  int seed;         // initial random seed

  /* --- hyperparameters --- */

  vec epsilon_range; 
  double nu_2;        // nu_2 = p + 2
  double nu_1;   // nu_1 = p + 2
  // Hyperparameter of the Inverse Wishart on Psi_1
  mat Psi_2; //  = eye<mat>(p,p);  
  // mean of the Normal prior on m_1
  vec m_2; // (p);  m_2.zeros();    
  // Covariance of the Normal prior on m_1
  mat inv_S_2; // =  eye<mat>(p,p); S_2 = S_2/1000;   
  // k_0 prior parameters
  vec tau_k0;  //  tau.fill(4);
  // alpha parameters
  vec tau_alpha;  // (2); tau_alpha.fill(1);
  // rho parameters
  vec rho_pm;  //(2); rho_0.fill(0.0);
  vec tau_rho;  // (2)
  // varphi parameters
  vec tau_varphi; // (2)
  vec varphi_pm; // (2); varphi_0.fill(0);  
  bool merge_step;  // 
  double merge_par;
  // latent indicator initial values
  uvec Z_input; 
  
  int length_chain; 
  
  vec saveRho, saveK0, saveEpsilon, saveVarphi;
  mat saveAlpha, saveW0;
  cube saveW1, saveMu, saveMu0, saveOmega;
  umat saveZ, saveS;
  
  
  void main_loop();
  Rcpp::List InitMuSigma(uvec Z, int k);

  Rcpp::List GenerateZetas( mat loglike,
                            arma::vec logW_0,
                            arma::mat logW_1,
                            double rho,
                            int n_cores );

  Rcpp::List UpdateZetas(   arma::cube mu, 
                            arma::cube Omega, 
                            arma::vec logW_0,
                            arma::mat logW_1,
                            double rho,
                            int n_cores=1 );
  
  arma::vec swapStep(   arma::mat N,
                        arma::vec alpha  );                          

  arma::vec UpdateAlpha(arma::vec alpha, arma::mat N, arma::vec alpha_par);
  
  double UpdateRho( arma::mat N );
  
  Rcpp::List UpdateLogWs(   arma::mat N, 
                            arma::vec alpha  );
                            
  Rcpp::List UpdateSMuSigma(  arma::uvec Z,
                              int k,  
                              arma::vec mu_0,
                              double varphi,
                              arma::mat Sigma_1, 
                              arma::mat Omega_1, 
                              double k_0, 
                              double epsilon,
                              arma::vec m_1 );      
                              
  double UpdateK0(  arma::cube Omega, 
                    arma::mat mu_0,
                    arma::vec m_1  );        
                    
  arma::mat UpdateSigma1(arma::cube Omega);
  
  arma::vec UpdateM1(  double k_0, 
                        arma::mat mu_0, 
                        arma::cube Omega );
                        
  double UpdateEpsilon( double epsilon,
                        arma::uvec S, 
                        arma::mat mu_0, 
                        arma::cube mu, 
                        arma::cube Omega,
                        double epsilon_par  );                      
                  
  double UpdateVarphi( arma::uvec S );
  
  Rcpp::List PriorSMuSigma(   double varphi,
                            arma::mat Sigma_1, 
                            arma::mat Omega_1, 
                            double k_0, 
                            double epsilon,
                            arma::vec m_1 ); 

  
  public:
  
  // constructor 
  MCMC( mat Y, 
        vec C, 
        Rcpp::List prior,
        Rcpp::List mcmc,
        Rcpp::List state );
        
  Rcpp::List get_chain();
      
  
};



#endif


