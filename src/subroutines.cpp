#include "approximations.h"

/* update_projection_sparse: update posterior of B for X.
 * num_factors: rank K
 * X: omics features that regress on the constructed factors.
 * Xobs: indicator matrix of observing value or not
 * weights: column-wise weights
 * nu_mat:  expected inverse variance of the noise.
 * meanFactors: expectation of the N by K factors
 * U2: second moments of the N by K factors
 * post_tmu: current posterior mean of non-zero B.
 * post_tsigma2: current posterior variance of non-zero B.
 * post_tpi: current posterior probability of non-zero B.
 * tau: prior variance of inactive B.
 * log_pi: expectation of prior log selection probability (for a path group).
 * log_minus_pi: expectation of prior log (1 - selection probability) (for a path group).
 * TODO: update weights_case
 * Done.
 */
double update_projection_sparse(const int num_factors, arma::mat& X,   const arma::mat& Xobs,
                                const arma::vec weights, const arma::vec weights_case,
                                arma::mat& nu_mat,
                                arma::mat& meanFactors,  arma::mat& U2,
                                arma::mat& post_tmu, arma::mat& post_tsigma2,  arma::mat& post_tpi,
                                arma::mat& tau, arma::mat& log_pi,  arma::mat& log_minus_pi){
  int p = X.n_cols;
  double Delta = 0.0;
  for(int j = 0; j < p; j++){
    arma::uvec obs = arma::find(Xobs.col(j) == 1);
    arma::vec x =  X.col(j);
    x = x(obs);
    arma::vec weight_obs = weights_case(obs);
    arma::vec nu_vec = nu_mat.col(j);
    nu_vec = nu_vec(obs);
    int n_obs = nu_vec.n_elem;
    arma::mat meanFactors_sub = meanFactors.rows(obs);
    arma::mat U2_sub = U2.rows(obs);
    arma::vec response = arma::zeros<arma::vec>(n_obs);
    response = response + x;
    /*
     * For computational efficiency, we keep track of the residuals.
     * When updating for factor k, the response is 
     *    residual' = residual+ current contribution from factor k.
     * After updating factor k, we also update the residual as
     *    residual = residual' - current contribution from factor k. 
     */
    for(int k = 0; k < num_factors; k++){
      double bkl = post_tmu(j,k) * post_tpi(j,k);
      response = response - meanFactors_sub.col(k) *bkl;
    }
    for(int k = 0; k < num_factors; k++){
      double bkl = post_tmu(j,k) * post_tpi(j,k);
      double bkl2 = post_tmu(j,k) * post_tmu(j,k)+post_tsigma2(j,k);
      arma::vec f = meanFactors_sub.col(k);
      arma::vec u = U2_sub.col(k);
      response = response + f * bkl;
      //update post_tmu, post_tsigma2 and post_tpi
      arma::vec response1 = (response % nu_vec % weight_obs) *weights(j);
      //arma::vec response1 = (response % nu_vec % weight_obs);
      //denominator of eq. (1)
      double tmp = arma::sum(nu_vec % u % weight_obs) *weights(j)+  tau(j,k);
      //double tmp = arma::sum(nu_vec % u % weight_obs)+  tau(j,k);
      //numerator of eq. (1)
      double tmp1 = sum(response1 % f);
      //deduct current contribution from B_{jk} to the EBLO using eq (1).
      //double delta1 = post_tpi(j,k) * (weights(j)*log_pi(j,k) -log(post_tpi(j,k)));
      double delta1 = post_tpi(j,k) * (log_pi(j,k) -log(post_tpi(j,k)));
      //double delta2 = (1.0 - post_tpi(j,k)) * (weights(j)*log_minus_pi(j,k) -log(1.0 - post_tpi(j,k)) - 0.5);
      double delta2 = (1.0 - post_tpi(j,k)) * (log_minus_pi(j,k) -log(1.0 - post_tpi(j,k)) - 0.5);
      //double delta2 = (1.0 - post_tpi(j,k)) * (log_minus_pi(j,k) -log(1.0 - post_tpi(j,k)));
      if(post_tpi(j,k) <= 1e-8){
        delta1 = 0.0;
      }
      if(post_tpi(j,k) >= (1.0-1e-8)){
        delta2 = 0.0;
      }
      Delta = Delta - (0.5 * post_tpi(j,k) *(-tmp * bkl2 + 2* tmp1 * post_tmu(j,k) + log(post_tsigma2(j,k) * tau(j,k))) +  delta1+ delta2);
      //update post_tmu, post_tsigma2, post_tpi with formulas (1)-(3)
      post_tsigma2(j,k)  = 1.0/tmp;
      post_tmu(j,k) = tmp1 * post_tsigma2(j,k);
      // double tmp2 = 0.5 *(post_tmu(j,k)* post_tmu(j,k)) / post_tsigma2(j,k) +
      //   0.5 * std::log(tau(j,k) * post_tsigma2(j,k)) +  weights(j)* log_pi(j,k) - weights(j)*log_minus_pi(j,k);
      double tmp2 = 0.5 *(post_tmu(j,k)* post_tmu(j,k)) / post_tsigma2(j,k) +
        0.5 * std::log(tau(j,k) * post_tsigma2(j,k)) +  log_pi(j,k) -log_minus_pi(j,k);
      if(tmp2 < 0){
        post_tpi(j,k) = std::exp(tmp2);
        post_tpi(j,k) = post_tpi(j,k)/(1.0+post_tpi(j,k));
      }else{
        post_tpi(j,k) = std::exp(-tmp2);
        post_tpi(j,k) = 1.0/(1.0+post_tpi(j,k));
      }
      //calculate the new ELBO (partial)
      response = response - f * bkl;
      bkl = post_tmu(j,k) * post_tpi(j,k);
      bkl2 = post_tmu(j,k) * post_tmu(j,k)+post_tsigma2(j,k);
      //delta1 = post_tpi(j,k) * (weights(j)*log_pi(j,k) -log(post_tpi(j,k)));
      //delta2 = (1.0 - post_tpi(j,k)) * (weights(j)*log_minus_pi(j,k) -log(1.0 - post_tpi(j,k)) - 0.5);
      delta1 = post_tpi(j,k) * (log_pi(j,k) -log(post_tpi(j,k)));
      delta2 = (1.0 - post_tpi(j,k)) * (log_minus_pi(j,k) -log(1.0 - post_tpi(j,k)) - 0.5);
      //delta2 = (1.0 - post_tpi(j,k)) * (log_minus_pi(j,k) -log(1.0 - post_tpi(j,k)));
      if(post_tpi(j,k) <= 1e-8){
        delta1 = 0.0;
      }
      if(post_tpi(j,k) >= (1.0-1e-8)){
        delta2 = 0.0;
      }
      Delta = Delta + (0.5 * post_tpi(j,k) *(-tmp * bkl2 + 2* tmp1 * post_tmu(j,k) + log(post_tsigma2(j,k) * tau(j,k)))
                         +  delta1+ delta2);    
      }
    
  }
  Delta = (Delta/(num_factors*p*1.0));
  return Delta;
}


double update_marginal_prob(arma::vec& x, arma::vec& nu_vec, arma::vec& meanfac, arma::vec& u2,
                            double taux, double log_pix, double log_minus_pix){
  double tmp = arma::sum(nu_vec % u2) +  taux;
  double tmp1 = sum(x % meanfac);
  double post_tsigma2x = 1.0/tmp;
  double post_tmux = tmp1 * post_tsigma2x;
  //update post_tmu, post_tsigma2, post_tpi with formulas (1)-(3)
  double post_tpix = 0.0;
  double tmp2 = 0.5 *(post_tmux * post_tmux) / post_tsigma2x +
    0.5 * std::log( taux * post_tsigma2x) +  log_pix - log_minus_pix;
  if(tmp2 < 0){
    post_tpix = std::exp(tmp2);
    post_tpix = post_tpix/(1.0+post_tpix);
  }else{
    post_tpix = std::exp(-tmp2);
    post_tpix = 1.0/(1.0+post_tpix);
  }
  return post_tpix;
}

/*
 * update_projection_constrained: update posterior of B for X.
 * Contraint is incoporated via proximal gradient descent
 * num_factor: rank K.
 * Y: response matrix.
 * Xobs: indicator matrix of observing value or not
 * nuYmat:  expected inverse variance of the noise.
 * meanFactors: expectation of the N by K factors
 * U2: second moments of the N by K factors
 * post_tmu: current posterior mean of non-zero B.
 * post_tsigma2: current posterior variance of non-zero B.
 * tau: prior variance of inactive B.
 * lower: lower bound on posterior magnitudes.
 * TODO: update weights_case, weights0
 * Done.
 */

double update_projection_constrained(const int num_factors, arma::mat& Y, const arma::mat& Yobs, const arma::vec weights0,
                                     const arma::vec weights_case,
                                     arma::mat& nuYmat, arma::mat& meanFactors,  arma::mat& U2,
                                     arma::mat& post_tmu, arma::mat& post_tsigma2,
                                     arma::mat& tau, double lower){
  int p = Y.n_cols;
  double Delta = 0.0;
  for(int j = 0; j < p; j++){
    arma::vec y =  Y.col(j);
    arma::uvec obs = arma::find(Yobs.col(j) == 1);
    arma::vec nu_vec = nuYmat.col(j);
    arma::vec weight_obs = weights_case(obs);
    nu_vec = nu_vec(obs);
    y = y(obs);
    int n_obs = nu_vec.n_elem;
    arma::mat meanFactors_sub = meanFactors.rows(obs);
    arma::mat U2_sub = U2.rows(obs);
    arma::vec response = arma::zeros<arma::vec>(n_obs);
    response = response + y;
    //Calculate matrix A and post_sigma2.
    arma::mat cov = arma::zeros<arma::mat>(num_factors, num_factors);
    for(int k =0; k < num_factors; k++){
      arma::vec f = meanFactors_sub.col(k);
      for(int k1 = 0; k1 < num_factors; k1++){
        arma::vec f1 = meanFactors_sub.col(k1);
        if(k!=k1){
          cov(k, k1) = arma::sum(f % f1 % nu_vec % weight_obs)*weights0(j);
        }
      }
    }
    for(int k = 0; k < num_factors; k++){
      arma::vec u = U2_sub.col(k);
      cov(k,k) = arma::sum( u % nu_vec%weight_obs)*weights0(j)+tau(j,k) * weights0(j);
      post_tsigma2(j,k) = 1.0/(cov(k,k));
    }
    arma::mat V1, U1;
    arma::vec S1;
    arma::svd(U1, S1, V1, cov, "std");
    if(U1.n_elem == 0){
      stop("svd failed");
    }
    arma::mat A = inv(cov);
    //respose = B in the updating algorithm
    response =meanFactors_sub.t() * (response % nu_vec%weight_obs)*weights0(j);
    arma::vec mu0 = A * response;
    double l2norm = sqrt(arma::sum(mu0 % mu0));
    if (l2norm < lower){
      //proximal gradient
      mu0 = post_tmu.row(j).t() * 1.0;
      double step_size = 1.0/(2*max(S1));
      int max_iter = 10000;
      double tol = 1e-8;
      int it  = 0;
      double err = tol * 2;
      while(it < max_iter & err > tol){
        arma::vec grad = cov * mu0- response;
        arma::vec mu1 = mu0 - step_size * grad;
        l2norm = sqrt(arma::sum(mu1 % mu1));
        if(l2norm < lower){
          mu1 = mu1 * lower/l2norm;
        }
        err = sqrt(arma::sum(arma::square(mu1 - mu0)));
        it = it + 1;
        mu0 = mu1 * 1.0;
      }
    }
    post_tmu.row(j) = mu0.t() * 1.0;
  }
  return Delta;
}





/*
 * update_factor_one: update factors for a given pattern
 * TODO: update weights_case, weights0
 * Done.
 */
double update_factor_one(const int num_factors, arma::mat& Y,  const arma::mat& Yobs,
                         arma::mat& X,  const arma::mat& Xobs,
                         const arma::mat& Z, const arma::vec weights, const arma::vec weights0,
                         const arma::vec weights_case,
                         arma::mat& post_mu, arma::mat& post_sigma2, arma::mat& post_pi,
                         arma::mat& post_tmuX, arma::mat& post_tsigma2X,
                         arma::mat& post_tpiX, arma::mat& post_tmuY, arma::mat& post_tsigma2Y,
                         arma::mat& tau,  arma::mat& log_pi, arma::mat& log_minus_pi,
                         arma::mat& meanFactors, arma::mat& nuXmat,
                         arma::mat& nuYmat, arma::mat& updatingOrders){
  int py = Y.n_cols;
  int n = Y.n_rows;
  int px = X.n_cols;
  int pz = Z.n_cols;
  double Delta = 0.0;
  for(int k = 0; k < num_factors; k++){
    arma::vec response =  arma::zeros<arma::vec>(n);
    arma::vec s_vec = arma::zeros<arma::vec>(n);
    //Calculate Rjk=response and w_{j'}nu_{j'}E[B_{j'k}^2]= s_vec: loop over X
    for(int j = 0; j < px; j++){
      double bkl  = 0.0;
      double bkl2 = 0.0;
      arma::vec x = X.col(j);
      arma::vec nu_vec = nuXmat.col(j);
      arma::vec xobs = Xobs.col(j);
      arma::uvec obs = arma::find(xobs == 1);
      bkl = post_tmuX(j,k) * post_tpiX(j,k);
      bkl2 = post_tmuX(j,k) * post_tmuX(j,k) * post_tpiX(j,k)+post_tsigma2X(j,k) * post_tpiX(j,k);
      response(obs) = response(obs)+ (x(obs) %nu_vec(obs)) * ( bkl* weights(j));
      s_vec(obs) += weights(j) * bkl2 * nu_vec(obs);
      for(int k1 = 0; k1 < num_factors; k1++){
        arma::vec f = meanFactors.col(k1);
        if(k1 != k){
          double cross = 0.0;
          cross = post_tmuX(j,k1) * post_tmuX(j,k) * post_tpiX(j,k1) *post_tpiX(j,k);
          response(obs) = response(obs) - (f(obs) % nu_vec(obs)) * (cross* weights(j));
        }
      }
    }
    //Calculate Rjk for factor k, feature j: loop over Y
    for(int j = 0; j < py; j++){
      double bkl  = 0.0;
      double bkl2 = 0.0;
      arma::vec y = Y.col(j);
      arma::vec nu_vec = nuYmat.col(j);
      arma::vec yobs = Yobs.col(j);
      arma::uvec obs = arma::find(yobs == 1);
      bkl = post_tmuY(j,k);
      bkl2 = post_tmuY(j,k) * post_tmuY(j,k)+post_tsigma2Y(j,k);
      response(obs) = response(obs)+ (y(obs) %nu_vec(obs)) *bkl* weights0(j);
      s_vec(obs) +=weights0(j)*bkl2 * nu_vec(obs);
      for(int k1 = 0; k1 < num_factors; k1++){
        arma::vec f = meanFactors.col(k1);
        if(k1 != k){
          double cross = 0.0;
          cross = post_tmuY(j,k1) * post_tmuY(j,k);
          response(obs) = response(obs) - (f(obs) % nu_vec(obs)) * cross*weights0(j);
        }
      }
    }
    //update
    for(int j0 = 0; j0 < pz; j0++){
      int j = updatingOrders(j0,k);
      double betajk = post_mu(j,k) * post_pi(j,k);
      double betajk2 = post_mu(j,k) * post_mu(j,k)+post_sigma2(j,k);
      meanFactors.col(k) = meanFactors.col(k) - Z.col(j) * betajk;
      arma::vec response1 = response - s_vec % meanFactors.col(k);
      //tmp = denominator of mu_jk.
      //double tmp =  arma::sum(s_vec % arma::square(Z.col(j))%weights_case)+tau(j,k)*weights0(0);
      //double tmp_tau = 1.0;
      double tmp_tau =  arma::sum(s_vec % arma::square(Z.col(j))%weights_case)+tau(j,k)*weights0(0);
      double tmp =  arma::sum(s_vec % arma::square(Z.col(j))%weights_case)+tmp_tau*weights0(0);
      //tmp1 = numerator of mu_jk.
      double tmp1 = sum(response1 % Z.col(j)%weights_case);
      double delta1 = post_pi(j,k) * (log_pi(j,k)*weights0(0) -log(post_pi(j,k)));
      double delta2 = (1.0 - post_pi(j,k)) * ( log_minus_pi(j,k)*weights0(0)  -log(1.0 - post_pi(j,k)) - 0.5);
      //double delta2 = (1.0 - post_pi(j,k)) * ( log_minus_pi(j,k)*weights0(0)  -log(1.0 - post_pi(j,k)));
      if(post_pi(j,k) <= 1e-8){
        delta1 = 0.0;
      }
      if(post_pi(j,k) >= (1.0-1e-8)){
        delta2 = 0.0;
      }
      Delta = Delta - (0.5 * post_pi(j,k) *(-tmp * betajk2 + 2* tmp1 * post_mu(j,k) + log(post_sigma2(j,k) *  tmp_tau))
                         +  delta1 + delta2);
      post_sigma2(j,k) = 1.0/tmp;
      post_mu(j,k) = tmp1 *  post_sigma2(j,k);
      double tmp2 = 0.5 * (post_mu(j,k)* post_mu(j,k))/ post_sigma2(j,k) +
        0.5 * std::log( tmp_tau * post_sigma2(j,k)) + (log_pi(j,k) - log_minus_pi(j,k))*weights0(0);
      if(tmp2 < 0){
        post_pi(j,k) = std::exp(tmp2);
        post_pi(j,k) = post_pi(j,k)/(1.0+post_pi(j,k));
      }else{
        post_pi(j,k) = std::exp(-tmp2);
        post_pi(j,k) = 1.0/(1.0+post_pi(j,k));
      }
      betajk2 = post_mu(j,k) * post_mu(j,k)+post_sigma2(j,k);
      betajk = post_mu(j,k) * post_pi(j,k);
      delta1 = post_pi(j,k) * (log_pi(j,k)*weights0(0) -log(post_pi(j,k)));
      delta2 = (1.0 - post_pi(j,k)) * ( log_minus_pi(j,k)*weights0(0) -log(1.0 - post_pi(j,k)) - 0.5);
      //delta2 = (1.0 - post_pi(j,k)) * ( log_minus_pi(j,k)*weights0(0) -log(1.0 - post_pi(j,k)));
      if(post_pi(j,k) <= 1e-8){
        delta1 = 0.0;
      }
      if(post_pi(j,k) >= (1.0-1e-8)){
        delta2 = 0.0;
      }
      Delta = Delta + (0.5 * post_pi(j,k) *(-tmp * betajk2 + 2* tmp1 * post_mu(j,k) + log(post_sigma2(j,k) *  tmp_tau))
                         +  delta1 + delta2);
      meanFactors.col(k) += Z.col(j) * betajk;
    }
  }
  Delta = (Delta/(num_factors*(pz*weights(0)+py*weights0(0))*1.0));
  return Delta;
}


/*
 * update_factor: update the factor compositions U = Zbeta
 * num_factors: rank K
 * Y, Yobs : response, observing indicator
 * X, Xobs : omics, observing indicator
 * Z: to contruct U.
 * weights: weigths for different U.
 * pattern_samples, pattern_features: sample_feature group cut
 */
double update_factor(const int num_factors, arma::mat& Y,  const arma::mat& Yobs,
                     arma::mat& X,  const arma::mat& Xobs,
                     const arma::mat& Z,  const arma::vec weights, const arma::vec weights0,
                     const arma::vec weights_case,
                     const List& pattern_samples, const List& pattern_features,
                     arma::cube& post_mu, arma::cube& post_sigma2, arma::cube& post_pi,
                     arma::mat& post_tmuX, arma::mat& post_tsigma2X,
                     arma::mat& post_tpiX, arma::mat& post_tmuY, arma::mat& post_tsigma2Y,
                     arma::mat& tauZ, arma::mat& log_pi, arma::mat& log_minus_pi,
                     arma::mat& nuXmat, arma::mat& nuYmat,
                     arma::mat& meanFactors,  arma::mat& U2,
                     arma::mat& updatingOrders){
  double Delta = 0.0;
  int num_pattern = pattern_samples.size();
  for(int k = 0; k < num_pattern; k++){
    //update for each pattern separately
    arma::uvec ii = pattern_samples(k);
    arma::uvec jj = pattern_features(k);
    ii = ii - 1;
    jj = jj -1;
    int p_sub = jj.n_elem;
    const arma::mat Zsub = Z.submat(ii, jj);
    arma::mat post_mu_sub = post_mu.slice(k);
    arma::mat post_sigma2_sub = post_sigma2.slice(k);
    arma::mat post_pi_sub = post_pi.slice(k);
    post_mu_sub = post_mu_sub.rows(jj);
    post_sigma2_sub = post_sigma2_sub.rows(jj);
    post_pi_sub = post_pi_sub.rows(jj);
    arma::mat Ysub = Y.rows(ii);
    const arma::mat Yobs_sub = Yobs.rows(ii);
    arma::mat Xsub = X.rows(ii);
    const arma::mat Xobs_sub = Xobs.rows(ii);
    arma::mat nuYmat_sub = nuYmat.rows(ii);
    arma::mat nuXmat_sub = nuXmat.rows(ii);
    arma::mat meanFactors_sub = meanFactors.rows(ii);
    double delta2 = update_factor_one(num_factors, Ysub, Yobs_sub, Xsub, Xobs_sub,
                                      Zsub, weights, weights0, weights_case,
                                      post_mu_sub, post_sigma2_sub,post_pi_sub ,
                                      post_tmuX, post_tsigma2X, post_tpiX,
                                      post_tmuY, post_tsigma2Y,
                                      tauZ, log_pi, log_minus_pi,
                                      meanFactors_sub, nuXmat, nuYmat, updatingOrders);
    arma::mat U2_sub = U2calculation(num_factors, Zsub, meanFactors_sub,
                                     post_mu_sub, post_sigma2_sub, post_pi_sub);
    if(delta2 > 0){
      Delta = Delta+delta2;
    }
    meanFactors.rows(ii) = meanFactors_sub;
    U2.rows(ii) = U2_sub;
    for(int j = 0; j < p_sub; j++){
      int j1 = jj(j);
      for(int l = 0; l < num_factors; l++){
        post_mu(j1,l,k) = post_mu_sub(j,l);
        post_sigma2(j1,l, k) = post_sigma2_sub(j,l);
        post_pi(j1,l, k) = post_pi_sub(j,l);
      }
    }
  }
  return Delta;
}

/*
 * Update the prior distribution of coefficient inverse variance for inactive features.
 * TODO: update weights, weights0
 */
double tau_update(const int& num_factors, const arma::vec& weights, const arma::vec& weights0,
                  const List& pattern_features, const List& functional_path,
                  const double a0, const double b0,arma::cube& post_mu,
                  arma::cube& post_sigma2, arma::cube& post_pi,
                  arma::mat& post_tmu, arma::mat& post_tsigma2,arma::mat& post_tpi,
                  arma::mat& tauZ, arma::mat& post_a0, arma::mat& post_b0, const double L,
                  const double L2){
  //posterior distribution for tau
  int num_path = functional_path.size();
  int num_pattern = pattern_features.size();
  double Delta = 0.0;
  for(int l = 0; l < num_path; l++){
    arma::uvec path_index = functional_path(l);
    path_index -=1;
    int p0 = path_index.n_elem;
    for(int k = 0; k < num_factors; k++){
      double post_a = a0;
      double post_b = b0;
      //add up the effects from factors
      for(int q = 0; q < num_pattern; q++){
        arma::uvec jj = pattern_features(q);
        jj -= 1;
        arma::uvec path_q = intersect(jj, path_index);
        int p1 = path_q.n_elem;
        for(int j0 = 0; j0 < p1; j0++){
          int j = path_q(j0);
          //double tmp = weights0(0)*((post_mu(j, k, q) * post_mu(j, k, q) + post_sigma2(j,k,q))*post_pi(j,k,q) + (1-post_pi(j,k,q))*1.0/tauZ(j,k))* .5;
          double tmp = weights0(0)*((post_mu(j, k, q) * post_mu(j, k, q) + post_sigma2(j,k,q)))* .5;
          post_a = post_a +  weights0(0)* .5;
          post_b = post_b + tmp;
        }
      }
      //add up the effects from weights
      for(int j0 = 0; j0 < p0; j0++){
        int j = path_index(j0);
        post_a = post_a +    .5;
        //post_b = post_b +   ((post_tmu(j, k) * post_tmu(j, k) + post_tsigma2(j,k))*post_tpi(j,k)+(1-post_tpi(j,k))*1.0/tauZ(j,k)) * .5;
        post_b = post_b +   ((post_tmu(j, k) * post_tmu(j, k) + post_tsigma2(j,k))) * .5;
      }
      Delta = Delta - ((post_a-post_a0(l,k)) * (boost::math::digamma(post_a0(l,k))-log(post_b0(l,k)))
                         - (post_b-post_b0(l,k)) * post_a0(l,k)/post_b0(l,k) +
                           std::lgamma(post_a0(l,k)) - post_a0(l,k) * log(post_b0(l,k)));
      double L_est = post_a/post_b;
      double post_a_use = post_a;
      double post_b_use = post_b;
      if(L_est >= L){
       post_b_use = post_a_use/L;
      }
      if(L_est <= L2){
        post_b_use = post_a_use/L2;
      }
      for(int j0 = 0; j0 < p0; j0++){
        int j = path_index(j0);
        tauZ(j,k) = post_a_use/post_b_use;
      }
      post_a0(l,k) = post_a_use;
      post_b0(l,k) = post_b_use;
      Delta = Delta + ((post_a-post_a0(l,k)) * (boost::math::digamma(post_a0(l,k))-log(post_b0(l,k)))
                         - (post_b-post_b0(l,k)) * post_a0(l,k)/post_b0(l,k) +
                           std::lgamma(post_a0(l,k)) - post_a0(l,k) * log(post_b0(l,k)));
      
    }
  }
  //Rcout << tauZ(0,0) << '\n';
  return Delta;
}


/*
 * Update the prior distribution of sparsity levels.
 * TODO: update weights, weights0
 */
double pi_update(const int& num_factors, const arma::vec& weights, const arma::vec& weights0,
                 const List& pattern_features, const List& functional_path,
                 const double a1, const double b1,
                 arma::cube& post_mu, arma::cube& post_sigma2, arma::cube& post_pi,
                 arma::mat& post_tmu, arma::mat& post_tsigma2,arma::mat& post_tpi,
                 arma::mat& log_pi, arma::mat& log_minus_pi,arma::mat& post_a1, arma::mat& post_b1, const double alpha0){
  //posterior distribution for tau
  int num_path = functional_path.size();
  int num_pattern = pattern_features.size();
  double Delta =0.0;
  for(int l = 0; l < num_path; l++){
    arma::uvec path_index = functional_path(l);
    path_index -= 1;
    int p0 = path_index.n_elem;
    for(int k = 0; k < num_factors; k++){
      double weight_sum = 0.0;
      double total_sum = 0.0;
      //add up the effects from factors
      for(int q = 0; q < num_pattern; q++){
        arma::uvec jj = pattern_features(q);
        jj -= 1;
        arma::uvec path_q = intersect(jj, path_index);
        int p1 = path_q.n_elem;
        for(int j0 = 0; j0 < p1; j0++){
          int j = path_q(j0);
          total_sum = total_sum + weights0(0);
          weight_sum = weight_sum + post_pi(j,k,q)*weights0(0);
        }
      }
      //add up the effects from weights
      for(int j0 = 0; j0 < p0; j0++){
        int j = path_index(j0);
        total_sum = total_sum + 1.0;
        weight_sum = weight_sum +post_tpi(j,k);
      }
      double post_a = a1+weight_sum ;
      double post_b = total_sum+b1 - weight_sum;
      //Rcout <<alpha0<<","<< post_a_use<<"," << post_b_use << '\n';
      double tmp1 = boost::math::digamma(post_a1(l,k)) - boost::math::digamma(post_a1(l,k)+post_b1(l,k));
      double tmp2 = boost::math::digamma(post_b1(l,k)) - boost::math::digamma(post_a1(l,k)+post_b1(l,k));
      double tmp3 = std::lgamma(a1) + std::lgamma(b1) - std::lgamma(a1+b1);
      double tmp4 = std::lgamma(post_a1(l,k)) + std::lgamma(post_b1(l,k)) - std::lgamma(post_a1(l,k)+post_b1(l,k));
      Delta = Delta - ((post_a - post_a1(l,k)) * tmp1 + (post_b- post_b1(l,k)) * tmp2 -
        (post_a+post_b) * tmp3 + tmp4);
      post_a1(l,k) = a1+weight_sum ;
      post_b1(l,k) = total_sum+b1 - weight_sum;
      double alpha1 = alpha0 * 0.1;
      if(alpha1 > 0.01){
        alpha1 = 0.01;
      }
      if(post_a>alpha0 * (post_a+post_b)){
        post_a1(l,k) = alpha0 * (post_a+post_b);
        post_b1(l,k) = (1.0-alpha0) * (post_a+post_b);
      }else if(post_a < alpha1*(post_a+post_b)){
        post_a1(l,k) = alpha1 * (post_a+post_b);
        post_b1(l,k) = (1.0-alpha1) * (post_a+post_b);
      }
      //Rcout << alpha0 << ',' << post_a1(l,k) << ',' <<post_b1(l,k)<<'\n';
      tmp1 = boost::math::digamma(post_a1(l,k)) - boost::math::digamma(post_a1(l,k)+post_b1(l,k));
      tmp2 = boost::math::digamma(post_b1(l,k)) - boost::math::digamma(post_a1(l,k)+post_b1(l,k));
      tmp3 = std::lgamma(a1) + std::lgamma(b1) - std::lgamma(a1+b1);
      tmp4 = std::lgamma(post_a1(l,k)) + std::lgamma(post_b1(l,k)) - std::lgamma(post_a1(l,k)+post_b1(l,k));
      Delta = Delta + ((post_a - post_a1(l,k)) * tmp1 + (post_b- post_b1(l,k)) * tmp2 -
        (post_a+post_b) * tmp3 + tmp4);
      for(int j0 = 0; j0 < p0; j0++){
        int j = path_index(j0);
        log_pi(j,k) = tmp1;
        log_minus_pi(j,k) = tmp2;
      }
    }
  }
  return Delta;
}

/*
 * Update Noise variances.
 * TODO: update weights, weights_case.
 * Done.
 */
double nu_update(arma::mat& Y, const arma::mat& Yobs, const arma::vec& weights, 
                 const arma::vec& weights_case,
                 arma::mat& meanFactors, arma::mat& U2,arma::mat& nu_mat,
                 const double a2, const double b2,
                 arma::mat& post_tmu, arma::mat& post_tsigma2,arma::mat& post_tpi,
                 arma::vec& post_a2, arma::vec& post_b2, bool sparse){
  int n = Y.n_rows;
  int p = Y.n_cols;
  arma::mat product = post_tmu % post_tpi;
  for(int j = 0; j < p; j++){
    product.row(j) = post_tmu.row(j);
  }
  arma::mat R = Y - meanFactors * product.t();
  arma::vec ones =arma::ones<arma::vec>(n);
  int num_factors = meanFactors.n_cols;
  arma::mat U20 = arma::zeros<arma::mat>(n, num_factors);
  for(int k = 0; k < num_factors; k++){
    U20.col(k) = arma::square(meanFactors.col(k));
  }
  double Delta = 0.0;
  for(int j = 0; j < p; j++){
    double correction = 0.0;
    arma::uvec obs = arma::find(Yobs.col(j) == 1);
    arma::vec r = R.col(j);
    //calculate LambdaU, LambdaB
    arma::mat LambdaU = arma::zeros<arma::mat>(num_factors,num_factors);
    arma::mat LambdaB = arma::zeros<arma::mat>(num_factors,num_factors);
    for(int k = 0; k < num_factors; k++){
      double tmp1 = 0.0;
      double tmp2 = 0.0;
      arma::vec u2 = U2.col(k);
      arma::vec u20 = U20.col(k);
      LambdaU(k,k) = arma::sum(u2(obs)%weights_case(obs))-arma::sum(u20(obs)%weights_case(obs));
      if(sparse){
        tmp1 = post_tmu(j,k)  * post_tpi(j,k);
        tmp2 = (post_tmu(j,k) * post_tmu(j,k)+post_tsigma2(j,k)) * post_tpi(j,k);
      }else{
        tmp1 = post_tmu(j,k);
        tmp2 = post_tmu(j,k) * post_tmu(j,k)+post_tsigma2(j,k);
      }
      LambdaB(k,k) = tmp2 - tmp1*tmp1;
      if(!sparse){
        //LambdaU(k,k)= 0.0;
      }
      correction += LambdaB(k,k) * arma::sum(u20(obs)%weights_case(obs));
      correction += tmp2 * LambdaU(k,k);
      //correction = correction + arma::sum(u2(obs)) * tmp2 - tmp1*tmp1*arma::sum(u20(obs));
    }
    if(!sparse){
      //correction = 0.0;
    }
    double tmp =  arma::sum(arma::square(r(obs))%weights_case(obs)) + correction;
    double post_a = a2 + weights(j) * arma::sum(weights_case(obs))/2.0;
    double post_b = b2 + weights(j) *tmp/2.0;
    Delta = Delta - ((post_a-post_a2(j)) * (boost::math::digamma(post_a2(j))-log(post_b2(j)))
                       - (post_b-post_b2(j)) * post_a2(j)/post_b2(j) +
                         std::lgamma(post_a2(j)) - post_a2(j) * log(post_b2(j)));
    nu_mat.col(j) = post_a/post_b * ones;
    double check = arma::min(nu_mat.col(j));
    if(check < 1e-20){
      stop("inadmissable variance estimation");
    }
    post_a2(j) = post_a;
    post_b2(j) = post_b;
    Delta = Delta +  (std::lgamma(post_a2(j)) - post_a2(j) * log(post_b2(j)));
  }
  return Delta;
}


