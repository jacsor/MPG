#include <RcppArmadillo.h>
#include "helpers.h"   
#include "gibbs.h"  

using namespace Rcpp;
using namespace arma;
using namespace std;


MCMC::MCMC( mat Y, 
            vec C, 
            Rcpp::List prior,
            Rcpp::List mcmc,
            Rcpp::List state ) :
            Y(Y), 
            C(C)
{
       
      p = Y.n_cols;  
      n = Y.n_rows;  
      J = C.max() + 1;
      K_0 = Rcpp::as<int>(prior["K_0"]);
      K_1 = Rcpp::as<int>(prior["K_1"]);
      num_iter = Rcpp::as<int>(mcmc["nskip"]) * Rcpp::as<int>(mcmc["nsave"]);
      num_burnin = Rcpp::as<int>(mcmc["nburn"]);
      num_thin = Rcpp::as<int>(mcmc["nskip"]);    
      num_display = Rcpp::as<int>(mcmc["ndisplay"]);    
      
      epsilon_range = Rcpp::as<vec>(prior["epsilon_range"]);
      m_2 = Rcpp::as<vec>(prior["m_2"]);
      nu_2 = Rcpp::as<double>(prior["nu_2"]);    
      nu_1 = Rcpp::as<double>(prior["nu_1"]);    
      Psi_2 = Rcpp::as<mat>(prior["Psi_2"]);
      inv_S_2 = inv(Rcpp::as<mat>(prior["S_2"]));
      tau_k0 = Rcpp::as<vec>(prior["tau_k0"]);
      tau_alpha = Rcpp::as<vec>(prior["tau_alpha"]);
      tau_rho = Rcpp::as<vec>(prior["tau_rho"]);
      rho_pm = Rcpp::as<vec>(prior["point_masses_rho"]);
      tau_varphi = Rcpp::as<vec>(prior["tau_varphi"]);
      varphi_pm = Rcpp::as<vec>(prior["point_masses_varphi"]);
      merge_step = Rcpp::as<bool>(prior["merge_step"]);
      merge_par = Rcpp::as<double>(prior["merge_par"]);
      Z_input = Rcpp::as<uvec>(state["Z"]);

      length_chain =  num_iter/num_thin;
      saveVarphi.set_size(length_chain);
      saveRho.set_size(length_chain);
      saveAlpha.set_size(length_chain,J+1);
      saveW0.set_size(length_chain,K_0);
      saveW1.set_size(J,K_1,length_chain);
      saveK0.set_size(length_chain);
      saveEpsilon.set_size(length_chain);
      saveS.set_size(length_chain,K_0+K_1);
      saveZ.set_size(length_chain,n);
      saveMu.set_size(J,(K_0+K_1)*p,length_chain);
      saveMu0.set_size(p,K_0+K_1,length_chain);
      saveOmega.set_size(p,(K_0+K_1)*p,length_chain);     
      
      main_loop();
                                  
}





void MCMC::main_loop()
{    
  
  /* --- support variables --- */
  
  // counter
  int km = 0;
  // latent indicators
  uvec S(K_0 + K_1);
  // number of observations per group and per component
  mat N(J,K_0 + K_1);
  // used in the swap step
  mat temp;
  vec indices;
  // latent assignments
  uvec Z = Z_input;

  /* --- parameters --- */
  
  // probability of misalignment and shift    
  double varphi = tau_varphi(0)/sum(tau_varphi);
  // proportion of the common mixture weights
  double rho = tau_rho(0)/sum(tau_rho);
  // link between mean and covariance
  double k_0 = 1;
  // perturbation parameter for the mean
  double epsilon  = mean(epsilon_range); //   epsilon(0) = 0.001; epsilon(1) = 1.0;
  double epsilon_old = epsilon;
  double epsilon_par = sqrt(K_0 + K_1);
  int epsilon_count = 0; 
  int epsilon_tot = 100;
  // mass parameter for the dirichlet prior on the mixture weights
  vec alpha(J+1); alpha.fill(1);
  vec alpha_old = alpha;
  vec alpha_par(J+1); alpha_par.fill(sqrt(K_0 + K_1));
  vec alpha_count(J+1); alpha_count.fill(0); 
  int alpha_tot = 100; 
  // mixture weights
  vec logW_0(K_0); logW_0.fill( log(1.0/K_0));
  mat logW_1(J, K_1); logW_1.fill( log(1.0/K_1));
  // mean \mu_{j,k}
  cube mu(J,p, K_0 + K_1);
  // centering of mean locations \mu_k
  mat mu_0(p,K_0 + K_1);  
  // covariance locations
  cube Sigma(p,p,K_0 + K_1);
  // precision locations
  cube Omega(p,p,K_0 + K_1);
  // centering for the Wishart prior 
  mat Sigma_1 = 10.0*eye<mat>(p,p); 
  mat Omega_1 = inv(Sigma_1);
  // mean of the mean
  vec m_1 = mean(Y,0).t(); 
  
  /* --- chain initialization --- */
    
  for(int k = 0; k < K_0 + K_1; k++ )
  {
    List tempMuSigma = InitMuSigma(Z, k); 
    mu.slice(k) = Rcpp::as<mat>(tempMuSigma["mu"]); 
    mu_0.col(k) = Rcpp::as<vec>(tempMuSigma["mu_0"]);
    Omega.slice(k) = Rcpp::as<mat>(tempMuSigma["Omega"]);   
    Sigma.slice(k) = Rcpp::as<mat>(tempMuSigma["Sigma"]);   
  }
  
  // assign each observation to a mixture component
  
  List tempZ = UpdateZetas(mu, Omega, logW_0, logW_1, rho, 1);   
  Z = Rcpp::as<uvec>(tempZ["Z"]);  
  N = Rcpp::as<mat>(tempZ["N"]);  
  
      
  /* --- let the chain run --- */

  for(int it=0; it<(num_iter + num_burnin); it++)
  {          
    if((it+1)%num_display==0)
      cout << "Iteration: " << it + 1 << " of " << num_iter + num_burnin << endl;
      
      
    indices = swapStep( N, alpha );
    
    if(indices(0) != indices(1))
    { 
      mu_0.swap_cols(indices(0),indices(1));    
      N.swap_cols(indices(0),indices(1));      
      temp = Omega.slice(indices(0));
      Omega.slice(indices(0)) = Omega.slice(indices(1));
      Omega.slice(indices(1)) = temp;         
      temp = mu.slice(indices(0));
      mu.slice(indices(0)) = mu.slice(indices(1));
      mu.slice(indices(1)) = temp;                                           
    }    
    
    // MERGE STEP
    if(merge_step)
    {
      for( int k=0; k< K_0 + K_1 - 1 ; k++ )
      {
        if( sum(N.col(k)) > 0 )
        {
          for(int kk=k+1; kk < K_0 + K_1 ; kk++)
          {
            if( sum(N.col(kk)) > 0  )
            {
              double kl_div = KL( mu_0.col(k), 
                                  mu_0.col(kk), 
                                  Sigma.slice(k), 
                                  Omega.slice(k), 
                                  Sigma.slice(kk), 
                                  Omega.slice(kk) ) / epsilon ;
              if( kl_div < R::qchisq(merge_par, (double)p, 1, 0) )
              {
                N.col(k) = N.col(k) + N.col(kk);
                N.col(kk) = zeros<vec>(J);
                
                List tempSMuSigma = PriorSMuSigma(   varphi,
                                                     Sigma_1, 
                                                     Omega_1, 
                                                     k_0, 
                                                     epsilon,
                                                     m_1 );                                                    
                S(kk) = Rcpp::as<unsigned>(tempSMuSigma["S"]);
                mu.slice(kk) = Rcpp::as<mat>(tempSMuSigma["mu"]); 
                mu_0.col(kk) = Rcpp::as<vec>(tempSMuSigma["mu_0"]);
                Omega.slice(kk) = Rcpp::as<mat>(tempSMuSigma["Omega"]);   
                Sigma.slice(kk) = Rcpp::as<mat>(tempSMuSigma["Sigma"]);                                        
                        
              }
  
            }
  
          }
        }
  
      }
    }
    
    alpha_old = alpha;  
    alpha = UpdateAlpha( alpha, N, alpha_par );
    if( it <= num_burnin )
    {
      for( int j=0; j<=J; j++)
      {
        if( alpha(j) != alpha_old(j) )
          alpha_count(j)++;
      }

        
      if( (it+1)  % alpha_tot == 0)
      {
        for( int j=0; j<=J; j++)
        {
          if( alpha_count(j) < 30 )
            alpha_par(j) *= 1.1;
          if( alpha_count(j) > 50 )  
            alpha_par(j) *= 0.9;
          alpha_count(j) = 0;
        }
      }      
    }  
    else
    {
      for( int j=0; j<=J; j++)
      {
        if( alpha(j) != alpha_old(j) )
          alpha_count(j)++;
      }
    }
        
    rho = UpdateRho( N );
    
    List tempW = UpdateLogWs( N, alpha );
    logW_0 = Rcpp::as<vec>(tempW["logW_0"]);
    logW_1 = Rcpp::as<mat>(tempW["logW_1"]);
    
    List tempZ = UpdateZetas(mu, Omega, logW_0, logW_1, rho, 1);   
    Z = Rcpp::as<uvec>(tempZ["Z"]);  
    N = Rcpp::as<mat>(tempZ["N"]);  
    
    for(int k=0; k < K_0 + K_1; k++)
    {      
      uvec Z_k = arma::find(Z==k);
      List tempSMuSigma = UpdateSMuSigma( Z,
                                          k,
                                          mu_0.col(k),
                                          varphi,
                                          Sigma_1, 
                                          Omega_1,
                                          k_0, 
                                          epsilon,
                                          m_1 ); 
      S(k) = Rcpp::as<unsigned>(tempSMuSigma["S"]);
      mu.slice(k) = Rcpp::as<mat>(tempSMuSigma["mu"]); 
      mu_0.col(k) = Rcpp::as<vec>(tempSMuSigma["mu_0"]);
      Omega.slice(k) = Rcpp::as<mat>(tempSMuSigma["Omega"]);   
      Sigma.slice(k) = Rcpp::as<mat>(tempSMuSigma["Sigma"]);   
      
    }  
    
    k_0 =  UpdateK0(Omega, mu_0, m_1);
    
    Sigma_1 = UpdateSigma1(Omega);
    Omega_1 = inv_sympd(Sigma_1);
    
    m_1 = UpdateM1( k_0, mu_0, Omega );
    
    epsilon_old = epsilon;
    epsilon = UpdateEpsilon(epsilon, S, mu_0, mu, Omega, epsilon_par);
    if( it <= num_burnin )
    {
      if( epsilon != epsilon_old )
        epsilon_count++;
        
      if( (it+1)  % epsilon_tot == 0)
      {
        if( epsilon_count < 30 )
          epsilon_par *= 1.1;
        if( epsilon_count > 50 )  
          epsilon_par *= 0.9;
        epsilon_count = 0;
      }      
    }  
    else
    {
      if( epsilon != epsilon_old )
        epsilon_count++;
    }
  
    varphi = UpdateVarphi( S );
    
    if( (it+1 > num_burnin) && ((it+1) % num_thin == 0))
    {  
      // save chain
      saveVarphi(km) = varphi;
      saveRho(km) = rho;
      saveK0(km) = k_0;
      saveEpsilon(km) = epsilon;
      saveAlpha.row(km) = alpha.t();
      saveW0.row(km) = exp(logW_0).t();   
      saveW1.slice(km) = exp(logW_1);   
      saveOmega.slice(km) = reshape( mat(Omega.memptr(), Omega.n_elem, 1, false), p, (K_0 + K_1)*p);  
      saveMu.slice(km) = reshape( mat(mu.memptr(), mu.n_elem, 1, false), J, (K_0 + K_1)*p);   
      saveMu0.slice(km) = mu_0;
      saveS.row(km) = S.t();          
      saveZ.row(km) = Z.t();  
      
      km++;        
    }
      
      
  }
  
  cout << endl << "MH acceptance rate " << endl;
  cout << "epsilon: " << (double)epsilon_count/num_iter << endl;
  cout << "alpha: " << alpha_count.t() / (double)num_iter << endl << endl;
  
}     




Rcpp::List MCMC::InitMuSigma(uvec Z, int k)
{
  mat EmpCov = eye<mat>(p,p) + cov(Y);
  mat tempInit = mvrnormArma(1, zeros<vec>(p) , eye<mat>(p,p));
  mat mu_k(J,p);
  vec mu_0k(p);  
  mat Sigma_k(p,p);
  mat Omega_k(p,p);
  
  uvec Zk = arma::find(Z==k);
  if( Zk.n_elem > 0 )
  {
    mu_0k = mean(Y.rows(Zk),0).t(); 
    
    if( Zk.n_elem > (uint)p + 2 )
      Sigma_k = cov(Y.rows(Zk));  
    else
      Sigma_k = EmpCov;  
      
    Omega_k = inv(Sigma_k);
    for(int j=0; j<J; j++)  
        mu_k.row(j) = mean(Y.rows(Zk),0);           
  }
  else
  {
    mu_0k = mean(Y,0).t() + tempInit.row(k).t();
    Sigma_k = EmpCov;  
    Omega_k = inv(Sigma_k);

    if( k < K_0 )
    {
      for(int j=0; j<J; j++)  
        mu_k.row(j) = mean(Y,0).t();        
    }
    else
    {
      for(int j=0; j<J; j++)  
        mu_k.row(j) = mean(Y,0) + tempInit.row(k);            
    }
  }                 
  
  return Rcpp::List::create(  
  Rcpp::Named( "mu" ) = mu_k,
  Rcpp::Named( "mu_0" ) = mu_0k,
  Rcpp::Named( "Sigma" ) = Sigma_k,
  Rcpp::Named( "Omega" ) = Omega_k) ;   
  
}





Rcpp::List MCMC::UpdateZetas(   arma::cube mu, 
                                arma::cube Omega, 
                                arma::vec logW_0,
                                arma::mat logW_1,
                                double rho,
                                int n_cores )
{
  mat N(J, K_0 + K_1); N.fill(0);
  uvec  Z(n);  Z.fill(0);  
  double tot_log_like = 0.0;
  mat loglike(n, K_0 + K_1);
  NumericVector U = runif(n);  
  uvec index(1); 
  
  if( rho > 0)
  {
    for(int k = 0; k < K_0; k++)
    {    
      index(0) = k;
      for(int j=0; j < J; j++)
      {
        uvec C_j = arma::find(C==j);
          loglike.submat(C_j,  index) = dmvnrm_arma_precision(Y.rows(C_j), mu.slice(k).row(j), 
            Omega.slice(k)) + log(rho);      
      }  
    }  
  }
  if( rho < 1 )
  {
    for(int k = K_0; k < K_0 + K_1; k++)
    {    
      index(0) = k;
      for(int j=0; j < J; j++)
      {
        uvec C_j = arma::find(C==j);
          loglike.submat(C_j,  index) = dmvnrm_arma_precision(Y.rows(C_j), mu.slice(k).row(j), 
            Omega.slice(k)) + log(1 - rho);      
      }  
    }    
  }  
  
  omp_set_num_threads(n_cores);
  // #pragma omp parallel for schedule(static) 
  for(int i=0; i < n; i++)
  {
    vec prob =  loglike.row(i).t();
    prob.rows(K_0, K_0 + K_1 - 1) = prob.rows(K_0, K_0 + K_1 - 1) + logW_1.row(C(i)).t();
    prob.rows(0,K_0-1) = prob.rows(0,K_0-1) + logW_0;
    prob = exp(prob);
    if( rho == 1 ) 
      prob.rows( 0, K_0-1 ).fill(0);
    else if(rho == 0)
      prob.rows( K_0, K_0 + K_1 - 1 ).fill(0);
    vec probsum = cumsum(prob);
    double x = U(i) * sum(prob);
    bool not_assigned = true;
    for (int k = 0; (k < K_0 + K_1) && not_assigned; k++)
    {  
      if(x <= probsum(k))
      {
        Z(i) = k;
        // S(i) = (k - k%(K_0 + K_1))/(K_0 + K_1);  // from 0 to 1
        not_assigned = false;
      }     
    }      
  }
    
  for(int i=0; i< n; i++)
  {      
    N(C(i),Z(i))++;  
    tot_log_like += loglike(i,Z(i));        
  }
  
  return Rcpp::List::create(  Rcpp::Named( "Z" ) = Z,
                              Rcpp::Named( "N" ) = N,
                              Rcpp::Named( "tot_log_like" ) = tot_log_like ) ;                                
                              
}



                               
arma::vec MCMC::UpdateAlpha(arma::vec alpha, arma::mat N, arma::vec alpha_par)
{
  vec output = alpha;    
  double log_ratio = 0;
  vec alpha_new(J+1);
  double temp;  
  vec counts(J+1);
  counts(0) = accu( N.cols(0,K_0 - 1) );
  counts.rows(1,J) = sum( N.cols(K_0, K_0 + K_1-1), 1 );

  for(int j=0; j<(J+1); j++)
  {
    if(counts(j)==0)
      output(j) = rgammaBayes(tau_alpha(0), tau_alpha(1));
    else
    {
      temp = rgammaBayes(  pow( alpha(j), 2 ) * alpha_par(j) , 
                        alpha(j) * alpha_par(j) );
      log_ratio = 0;                  
      
      if( j == 0)
      {
        log_ratio += marginalLikeDirichlet( sum(N.cols(0,K_0-1),0).t(), (temp/K_0)*ones<vec>(K_0)  );
        log_ratio -= marginalLikeDirichlet( sum(N.cols(0,K_0-1),0).t(), (alpha(0)/K_0)*ones<vec>(K_0)  );
      }
      else
      {
        log_ratio += marginalLikeDirichlet( N.cols(K_0, K_0 + K_1 -1).row(j-1).t(), (temp/K_1)*ones<vec>(K_1)  );
        log_ratio -= marginalLikeDirichlet( N.cols(K_0, K_0 + K_1 -1).row(j-1).t(), (alpha(j)/K_1)*ones<vec>(K_1)  );
      }

      log_ratio += R::dgamma(temp, tau_alpha(0), 1/tau_alpha(1), 1);
      log_ratio -= R::dgamma(alpha(j), tau_alpha(0), 1/tau_alpha(1), 1);
      
      log_ratio += R::dgamma(alpha(j), pow(temp,2)* alpha_par(j)  , 
                           1/temp/alpha_par(j)  , 1);                          
      log_ratio -= R::dgamma(temp, pow(alpha(j),2)* alpha_par(j) , 
                           1/alpha(j)/alpha_par(j), 1);      
      
      // cout << exp(log_ratio) << " " ;
      if( exp(log_ratio) > R::runif(0,1) )
          output(j) = temp;
    }    
  }
  // cout << endl;
  return output;

}



double MCMC::UpdateRho( arma::mat N )
{
  if( rho_pm(0)==1 )
    return 0.0;
  else if( rho_pm(1)==1 )
    return 1.0;
  else
  {
    vec counts(2); counts.fill(0);
    counts(0) = accu( N.cols(0,K_0-1) );
    counts(1) = accu( N.cols(K_0,K_0+K_1-1) );
    for( int s=0; s<2; s++)
    {    
      if( ( counts(s) == n ) && ( rho_pm(s) > 0 ) )
      {
         double log_den = log_exp_x_plus_exp_y( log(rho_pm(s)) , 
          log(1 - sum(rho_pm) ) + marginalLikeDirichlet(counts, tau_rho) );
    
        if( exp( log(rho_pm(s)) - log_den  ) <   R::runif(0,1) )
          return as<double>(rbeta(1, tau_rho(0) + counts(0), tau_rho(1) + counts(1)));
        else
          return (1-s);
      }
    }
    return as<double>(rbeta(1, tau_rho(0) + counts(0), tau_rho(1) + counts(1)));
  }

}




Rcpp::List MCMC::UpdateLogWs(   arma::mat N, 
                                arma::vec alpha  )
{
  vec logW_0(K_0);
  mat logW_1(J,K_1);  
  logW_0 = rDirichlet( sum(N.cols(0,K_0-1),0).t() +  alpha(0) * ones<vec>(K_0) / K_0 );

  for(int j=0; j<J; j++)
    logW_1.row(j) = rDirichlet( N.cols(K_0,K_0+K_1-1).row(j).t() +  alpha(j+1) * ones<vec>(K_1) / K_1  ).t();
  
  return Rcpp::List::create(  Rcpp::Named( "logW_0" ) = logW_0,
                              Rcpp::Named( "logW_1" ) = logW_1 ) ;      
}





Rcpp::List MCMC::UpdateSMuSigma(  arma::uvec Z,
                                  int k, 
                                  arma::vec mu_0,
                                  double varphi,
                                  arma::mat Sigma_1, 
                                  arma::mat Omega_1, 
                                  double k_0, 
                                  double epsilon,
                                  arma::vec m_1 ) 
{ 
  uvec Z_k = arma::find(Z==k);  
  mat data_group = Y.rows(Z_k);
  vec C_k = C(Z_k);
  int p = data_group.n_cols;
  int N_k = data_group.n_rows;   
  unsigned S_k;
  mat Omega(p,p);
  mat Sigma(p,p);
  mat mu(p,J);
  mat mu_0new(p,1);  
  double r = R::runif(0,1);
  
  if(N_k == 0)
  {
    Omega = rWishartArma(Omega_1, nu_1);
    Sigma = inv_sympd( Omega );       
    mu_0new = trans(mvrnormArma(1, m_1, Sigma/k_0));    
    
    if( r < varphi )
    {
      S_k = 0;
      mu = repmat(mu_0, 1, J);
    }
    else 
    {
      S_k = 1;
      for(int j=0; j<J; j++)
        mu.col(j) = trans( mvrnormArma(1, mu_0new,  Sigma*epsilon ));      
    }
  }
  else  // N_k > 0
  {
    double sign_0, sign_1;
    mat Psi_0(p,p), Psi_1(p,p);
    double extra_piece_var_1 = 0;
    vec m1_0(p), m1_1(p);
    double log_extra_piece = 0;
    vec log_det_Psi(2); log_det_Psi.fill(0);
    
    vec n_jk(J);
    mat mean_jk(p,J); mean_jk.fill(0);
    vec mean_k = mean(data_group,0).t();
    mat SS_jk(p,p), ss_jk_1(p,p);
    
    vec log_prob(2); 
    double prob;
    double log_prob_sum;
    
    if( ( varphi > 0  ) && ( varphi < 1) )
    {
      // marginal likelihood under model 0
      mat mean_k_rep = repmat( mean_k.t(), N_k, 1);
      mat SS_k = ( data_group - mean_k_rep ).t() * ( data_group - mean_k_rep );
      mat ss_k = N_k * ( mean_k - mu_0 ) * ( mean_k - mu_0 ).t();        
      Psi_0 = inv_sympd( Sigma_1 + SS_k + ss_k );    
      log_det(log_det_Psi(0), sign_0, Psi_0); 
      m1_0 = ( N_k * mean_k + k_0 * m_1 ) / (N_k + k_0);
      
      log_prob(0) = (nu_1 + N_k)/2.0 * log_det_Psi(0) + log( varphi );

      // marginal likelihood under model 1 
      extra_piece_var_1 = k_0;
      m1_1 = k_0 * m_1;     
      SS_jk.fill(0);
      ss_jk_1.fill(0);

      for(int j=0; j<J; j++)
      {
        uvec indices = find(C_k==j);
        n_jk(j) = indices.n_elem;     
        
        if (n_jk(j) > 0)
        {
            mean_jk.col(j) = mean(data_group.rows(indices),0).t();
            mat mean_jk_rep = repmat( trans(mean_jk.col(j)),(int)n_jk(j), 1);
            SS_jk = SS_jk + (data_group.rows(indices) - mean_jk_rep).t() * ( data_group.rows(indices) - mean_jk_rep );
            ss_jk_1 = ss_jk_1 + (mean_jk.col(j) - mu_0) * (mean_jk.col(j) - mu_0).t() / (epsilon + 1.0/n_jk(j));
            log_extra_piece +=  log( epsilon*n_jk(j) + 1.0 );
            extra_piece_var_1 +=  n_jk(j) / (epsilon * n_jk(j) + 1.0);
            m1_1 = m1_1 +  n_jk(j) / (epsilon * n_jk(j) + 1.0) * mean_jk.col(j);
           
        }
      }          
      Psi_1 = inv_sympd( Sigma_1 + SS_jk + ss_jk_1 );
      log_det(log_det_Psi(1), sign_1, Psi_1);
      
      log_prob(1) = (nu_1 + N_k)/2.0 * log_det_Psi(1) - p/2.0 * log_extra_piece + log( 1.0 - varphi );
      
      log_prob_sum = log_exp_x_plus_exp_y( log_prob(0), log_prob(1) );
      prob = exp( log_prob(0) - log_prob_sum  );
    }
    else
      prob = varphi;
    
    if( r < prob )
    {
      S_k = 0;
      Omega = rWishartArma(Psi_0, nu_1 + N_k);
      Sigma = inv_sympd( Omega ); 
      mu_0new = trans( mvrnormArma(1, m1_0, Sigma/(k_0+N_k)) );      
      mu = repmat(mu_0new, 1, J);
    }
    else 
    {
      S_k = 1;
      Omega = rWishartArma(Psi_1, nu_1 + N_k);
      Sigma = inv_sympd( Omega ); 
      m1_1 = m1_1 / extra_piece_var_1;
      mu_0new = trans( mvrnormArma(1, m1_1, Sigma/extra_piece_var_1));
      for(int j=0; j<J; j++)
      {
        if( n_jk(j) > 0 )
        {
          mu.col(j) = trans(mvrnormArma(1, (n_jk(j)*mean_jk.col(j) + 1.0/epsilon*mu_0new)/(n_jk(j) + 1.0/epsilon), 
            Sigma/(n_jk(j) + 1.0/epsilon)));
        }
        else
          mu.col(j) = trans( mvrnormArma(1, mu_0new,  Sigma*epsilon) );           
      }      
    }
  }
           
  return Rcpp::List::create(  
    Rcpp::Named( "S" ) = S_k,
    Rcpp::Named( "mu" ) = mu.t(),
    Rcpp::Named( "mu_0" ) = mu_0new,
    Rcpp::Named( "Sigma" ) = Sigma,
    Rcpp::Named( "Omega" ) = Omega) ;    
};





double MCMC::UpdateK0(  arma::cube Omega, 
                        arma::mat mu_0,
                        arma::vec m_1  )
{
  double tau_2_tot = tau_k0(1);
  double tau_1_tot = tau_k0(0) + p*(K_0 + K_1);
  for(int k=0; k < (K_0 + K_1); k++)
      tau_2_tot += as_scalar( (mu_0.col(k) - m_1).t() * Omega.slice(k) * (mu_0.col(k) - m_1));  

  return rgammaBayes(tau_1_tot/2, tau_2_tot/2);
};






arma::mat MCMC::UpdateSigma1(arma::cube Omega)
{
  mat psi_2_tot = Psi_2;
  for(int k=0; k< (K_0 + K_1); k++)
      psi_2_tot += Omega.slice(k);

  return( rWishartArma(inv_sympd(psi_2_tot), (K_0 + K_1)*nu_1 + nu_2) );
  
};


arma::vec MCMC::UpdateM1(   double k_0, 
                            arma::mat mu_0, 
                            arma::cube Omega )
{
  mat precision = inv_S_2;
  vec meanM = inv_S_2*m_2;
  for(int k=0; k< K_0 + K_1; k++)
  {
      precision += k_0*Omega.slice(k);
      meanM += k_0 * ( Omega.slice(k)*mu_0.col(k) );
  }
  mat variance = inv_sympd(precision);
  mat output = mvrnormArma(1, variance*meanM, variance);
  return( output.row(0).t() );
};



double MCMC::UpdateEpsilon( double epsilon,
                            arma::uvec S, 
                            arma::mat mu_0, 
                            arma::cube mu, 
                            arma::cube Omega,
                            double epsilon_par  )
{
  
  double counts;
  uvec index = find(S==1);
  counts = index.n_elem;
  
  double epsilon_new;
  double numerator = 0;
  double denominator = 0;
  double output = epsilon;
  double eps_diff = epsilon_range(1) - epsilon_range(0);
  double a,b;
  
  if(counts==0)
    output = R::runif( epsilon_range(0), epsilon_range(1) );   // sample from the prior
  else
  {
    a = exp(  log( epsilon_par ) + log( epsilon - epsilon_range(0) ) - log( eps_diff ) ) ;
    b = epsilon_par - a ;
    
    epsilon_new = (eps_diff)*as<double>(rbeta(1, a, b)) + epsilon_range(0);
    
    denominator = dGeneralizedBeta(epsilon_new, a, b , epsilon_range  );
    
    a = exp(  log( epsilon_par ) + log( epsilon_new - epsilon_range(0) ) - log( eps_diff ) );
    b = epsilon_par - a ;
    
    numerator = dGeneralizedBeta(epsilon, a, b , epsilon_range  );
    
    for(int k=0; k< K_0 + K_1; k++)
    {
      if( S(k) == 1 )
      {
        numerator += sum(dmvnrm_arma_precision(mu.slice(k),mu_0.col(k).t(), 
          Omega.slice(k)/epsilon_new));
        denominator += sum(dmvnrm_arma_precision(mu.slice(k),mu_0.col(k).t(), 
          Omega.slice(k)/epsilon ));
      }
    }      
    // cout << exp(numerator - denominator) << endl;
    
    if( exp(numerator - denominator) > R::runif(0,1) )
      output = epsilon_new;
    
  }
    
  return output;  
  
  
};







double MCMC::UpdateVarphi( arma::uvec S )
{
  if( varphi_pm(0)==1 )
    return 0.0;
  else if( varphi_pm(1)==1 )  
    return 1.0;
  else
  {
    vec counts(2); 
    uvec index = find(S==0);
    counts(0) = index.n_elem;
    index = find(S==1);
    counts(1) = index.n_elem;
        
    for(int s=0; s<2; s++)
    {      
      if( ( counts(s) == K_0 + K_1 ) && ( varphi_pm(s) > 0 ) )
      {
         double log_den = log_exp_x_plus_exp_y( log(varphi_pm(s)) , 
          log(1 - sum(varphi_pm)  ) + marginalLikeDirichlet(counts, tau_varphi) );
    
        if( exp( log(varphi_pm(s)) - log_den  ) <   R::runif(0,1) )
          return as<double>(rbeta(1, tau_varphi(0) + counts(0), tau_varphi(1) + counts(1)));
        else  
          return (1-s);         
      }    
    }
    return as<double>(rbeta(1, tau_varphi(0) + counts(0), tau_varphi(1) + counts(1)));
  }      
}



Rcpp::List MCMC::get_chain()
{
  return Rcpp::List::create(  
    Rcpp::Named( "alpha" ) = saveAlpha,
    Rcpp::Named( "epsilon" ) = saveEpsilon,
    Rcpp::Named( "k_0" ) = saveK0,
    Rcpp::Named( "mu" ) = saveMu,
    Rcpp::Named( "mu_0" ) = saveMu0,
    Rcpp::Named( "Omega" ) = saveOmega,
    Rcpp::Named( "rho" ) = saveRho,
    Rcpp::Named( "Z" ) = saveZ,
    Rcpp::Named( "S" ) = saveS,
    Rcpp::Named( "varphi" ) = saveVarphi,
    Rcpp::Named( "w0" ) = saveW0,
    Rcpp::Named( "w1" ) = saveW1
    );
};






arma::vec MCMC::swapStep(   arma::mat N,
                            arma::vec alpha  )
{
  vec probs = sqrt( sum(N,0).t()  );    
  vec unif_probs; 
  vec indices(2);   
  vec same_indices(2); same_indices.fill(0);
  
  indices(0) = sampling(probs);
  if(indices(0)  < K_0)
  {
    unif_probs = ones<vec>(K_1);
    indices(1) = sampling(unif_probs) + K_0;
  }
  else
  {
    unif_probs = ones<vec>(K_0);
    indices(1) = sampling(unif_probs);
  }    
           
  mat N_new = N;
  N_new.swap_cols(indices(0), indices(1));
  
  vec counts(2), counts_new(2);
  counts(0) = accu( N.cols(0,K_0-1) );
  counts(1) = accu( N.cols(K_0,K_0+K_1-1) );
  counts_new(0) = accu( N_new.cols(0,K_0-1) );
  counts_new(1) = accu( N_new.cols(K_0,K_0+K_1-1) );  
  
  double numerator = 0, denominator = 0;

  if( (counts_new(0) == 0) && (rho_pm(0)>0) )
    numerator += log_exp_x_plus_exp_y( log(rho_pm(0)), log(1 - sum(rho_pm)) + marginalLikeDirichlet(counts_new, tau_rho ));        
  else if( (counts_new(1) == 0) && (rho_pm(1)>0) )
    numerator += log_exp_x_plus_exp_y( log(rho_pm(1)), log(1 - sum(rho_pm)) + marginalLikeDirichlet(counts_new, tau_rho ));      
  else
    numerator += log(1 - sum(rho_pm)) + marginalLikeDirichlet( counts_new, tau_rho );      
    
  if( (counts(0) == 0) && (rho_pm(0)>0) )
    denominator += log_exp_x_plus_exp_y( log(rho_pm(0)), log(1 - sum(rho_pm)) + marginalLikeDirichlet(counts, tau_rho  ));      
  else if( (counts(1) == 0) && (rho_pm(1)>0) )
    denominator += log_exp_x_plus_exp_y( log(rho_pm(1)), log(1 - sum(rho_pm)) + marginalLikeDirichlet(counts, tau_rho  ));      
  else
    denominator += log(1 - sum(rho_pm)) + marginalLikeDirichlet( counts, tau_rho );      

  numerator += marginalLikeDirichlet( sum(N_new.cols(0,K_0-1),0).t(), alpha(0)*ones<vec>(K_0)/K_0 );
  denominator += marginalLikeDirichlet( sum(N.cols(0,K_0-1),0).t(), alpha(0)*ones<vec>(K_0)/K_0 );
  
  for(int j=0; j<J; j++)
  {
    numerator += marginalLikeDirichlet( N_new.cols(K_0,K_0+K_1-1).row(j).t(), alpha(j+1)*ones<vec>(K_1)/K_1 );
    denominator += marginalLikeDirichlet( N.cols(K_0,K_0+K_1-1).row(j).t(), alpha(j+1)*ones<vec>(K_1)/K_1 );
  }
  
  double acceptance_log_prob = numerator - denominator;
  
  if( exp(acceptance_log_prob) > R::runif(0,1) )
  {
    return indices;
  }
  else   
    return same_indices;
};






Rcpp::List MCMC::PriorSMuSigma(   double varphi,
                                  arma::mat Sigma_1, 
                                  arma::mat Omega_1, 
                                  double k_0, 
                                  double epsilon,
                                  arma::vec m_1 ) 
{ 
  unsigned S_k;
  mat Omega(p,p);
  mat Sigma(p,p);
  mat mu(p,J);
  mat mu_0new(p,1);  
  double r = R::runif(0,1);

  Omega = rWishartArma(Omega_1, nu_1);
  Sigma = inv_sympd( Omega );       
  mu_0new = trans(mvrnormArma(1, m_1, Sigma/k_0));    
    
  if( r < varphi )
  {
    S_k = 0;
    mu = repmat(mu_0new, 1, J);
  }
  else 
  {
    S_k = 1;
    for(int j=0; j<J; j++)
      mu.col(j) = trans( mvrnormArma(1, mu_0new,  Sigma*epsilon ));      
  }
  
  return Rcpp::List::create(  
    Rcpp::Named( "S" ) = S_k,
    Rcpp::Named( "mu" ) = mu.t(),
    Rcpp::Named( "mu_0" ) = mu_0new,
    Rcpp::Named( "Sigma" ) = Sigma,
    Rcpp::Named( "Omega" ) = Omega) ;    
};

