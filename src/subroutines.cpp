#include <RcppArmadillo.h>
#include <boost/math/special_functions/digamma.hpp>
#include <cmath>
// [[Rcpp::depends(BH)]]
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;



/*
 * Logistic binomial regression
 */
void UPXI_update_binomial(int idx, arma::mat& S, arma::mat& Y, arma::mat& Yapprox,
                          arma::vec& offset, arma::vec& Delta,
                          arma::vec& intercepts, arma::mat& UPXI, arma::mat& nu_mat,
                          const double robust_eps){
  //Jakolaa bound
  int n = UPXI.n_rows;
  int K = UPXI.n_cols;
  nu_mat.col(idx) = arma::zeros<arma::vec>(n);
  Yapprox.col(idx) = arma::zeros<arma::vec>(n);
  for(int k = 0; k < K; k++){
    arma::vec tmp = offset + intercepts(k);
    UPXI.col(k) = tmp % tmp + Delta;
    UPXI.col(k) = arma::sqrt(UPXI.col(k));
    arma::vec lambdas = exp(-UPXI.col(k));
    lambdas = (1.0 - lambdas)/(1.0+lambdas)/(4 * UPXI.col(k));
    for(int i = 0; i < n; i++){
      if(UPXI(i,k) < 1e-4){
        lambdas(i) = 1.0/8.0;
      }
    }
    lambdas = lambdas + robust_eps;
    nu_mat.col(idx) += 2.0 * lambdas;
    Yapprox.col(idx) += (2 * S.col(k) - 1.0)/4.0  - lambdas * intercepts(k);
  }
  Yapprox.col(idx) = 2.0/nu_mat.col(idx) % Yapprox.col(idx);
}

void update_binomial_approximation(arma::mat& Y, const arma::mat& Yobs,
                                   const arma::vec nclasses,  arma::mat& Yapprox,
                                   arma::mat& meanFactors, arma::mat& U2,
                                   arma::mat& post_tmu, arma::mat& post_tsigma2,
                                   arma::mat& post_tpi, List& intercepts, arma::mat& nu_mat,
                                   const double robust_eps){
  int p = Y.n_cols;
  int n = Y.n_rows;
  int num_factors = meanFactors.n_cols;
  for(int j = 0; j < p; j++){
    //number of classes
    int K = nclasses(j) - 1;
    arma::vec intercept_prev = intercepts(j);
    arma::mat S = arma::zeros<arma::mat>(n, K);
    arma::mat UPXI = arma::zeros<arma::mat>(n, K);
    for(int k = 0; k < K; k++){
      arma::vec tmp = arma::zeros<arma::vec>(n);
      for(int i = 0; i < n; i++){
        if(Y(i,j) >= (k+1)){
          S(i,k) = 1.0;
        }
      }
    }
    arma::vec b = (post_tmu.row(j) % post_tpi.row(j)).t();
    arma::vec b2 = ((post_tmu.row(j) % post_tmu.row(j) + post_tsigma2.row(j)) % post_tpi.row(j)).t();
    arma::vec offset = meanFactors * b;
    arma::vec Delta = arma::zeros<arma::vec>(n);
    for(int k = 0; k < num_factors; k++){
      Delta += U2.col(k) * b2(k) - (meanFactors.col(k) % meanFactors.col(k)) * b(k) * b(k);
    }
    for(int k = 0; k < K; k++){
      UPXI.col(k) = (offset + intercept_prev(k));
      UPXI.col(k)  = UPXI.col(k) % UPXI.col(k) + Delta;
      UPXI.col(k) = arma::sqrt(UPXI.col(k));
    }
    //update intercepts USING the non-missing values!
    arma::vec tmp1 = Yobs.col(j);
    for(int k = 0; k < K; k++){
      arma::vec lambdas = exp(-UPXI.col(k));
      lambdas = (1.0 - lambdas)/(1.0+lambdas)/(4 * UPXI.col(k))+robust_eps;
      intercept_prev(k) = arma::sum((2 * S.col(k) - 1.0) % tmp1)/4 - arma::sum(lambdas % offset % tmp1);
      intercept_prev(k) = intercept_prev(k)/arma::sum(lambdas % tmp1);
    }
    //update UPXI matrix and Y approx
    UPXI_update_binomial(j, S, Y,  Yapprox, offset,  Delta,
                         intercept_prev,  UPXI, nu_mat, robust_eps);
    intercepts(j) = intercept_prev;
  }
}





/*
 * Logistic ordinal regression
 */
void UPXI_update_ordinal(int idx,  arma::mat& Y,  arma::mat& Yapprox,
                         arma::vec& offset, arma::vec& Delta,
                         arma::vec& intercepts, arma::mat& UPXI, arma::mat& nu_mat,
                         const double robust_eps){
  //Jakolaa bound
  int n = UPXI.n_rows;
  int K = UPXI.n_cols;
  nu_mat.col(idx) = arma::zeros<arma::vec>(n);
  Yapprox.col(idx) = arma::zeros<arma::vec>(n);
  for(int k = 0; k < K; k++){
    arma::vec tmp = offset + intercepts(k);
    UPXI.col(k) = tmp % tmp + Delta;
    UPXI.col(k) = arma::sqrt(UPXI.col(k));
  }
  for(int i = 0; i < n; i++){
    if(Y(i,idx) == 0){
      double lambda = exp(-UPXI(i,0));
      if(UPXI(i,0) > 1e-4){
        lambda = (1.0 - lambda)/(1.0+lambda)/(4 * UPXI(i,0));
      }else{
        lambda = 1.0/8.0;
      }
      lambda = lambda + robust_eps;
      Yapprox(i, idx) = -intercepts(0) - 1.0/(4 * lambda);
      nu_mat(i, idx) = 2 * lambda;
    }else if(Y(i, idx) == K){
      double lambda = exp(-UPXI(i, K-1));
      if(UPXI(i,K-1) > 1e-4){
        lambda = (1.0 - lambda)/(1.0+lambda)/(4 * UPXI(i,0))+ robust_eps;
      }else{
        lambda = 1.0/8.0;
      }
      lambda = lambda + robust_eps;
      Yapprox(i, idx) = - intercepts(K-1) + 1.0/(4 * lambda);
      nu_mat(i, idx) = 2 * lambda;
    }else{
      double lambda1 = exp(-UPXI(i, Y(i,idx)-1));
      double lambda2 = exp(-UPXI(i, Y(i,idx)));
      if(UPXI(i,Y(i,idx)-1) > 1e-4){
        lambda1 = (1.0 - lambda1)/(1.0+lambda1)/(4 * UPXI(i,Y(i,idx)-1));
      }else{
        lambda1 = 1.0/8.0;
      }
      if(UPXI(i,Y(i,idx)) > 1e-4){
        lambda2 = (1.0 - lambda2)/(1.0+lambda2)/(4 * UPXI(i,Y(i,idx)));
      }else{
        lambda2 = 1.0/8.0;
      }
      lambda1 = lambda1 + robust_eps;
      lambda2 = lambda2 + robust_eps;
      nu_mat(i, idx) =  (lambda1 + lambda2);
      Yapprox(i, idx) = - intercepts(Y(i, idx)-1)* lambda1 - intercepts(Y(i, idx)) * lambda2;
      Yapprox(i, idx) = Yapprox(i, idx)/(nu_mat(i,idx));
      nu_mat(i, idx) = 2 * nu_mat(i, idx);
    }
  }
}



double obj_fun_ordinal(arma::vec& intercepts,
                       arma::vec& linear_coefs, arma::vec& quad_coefs, arma::vec& logDelta_coefs){
  int K = intercepts.n_elem;
  double obj_val = 0.0;
  arma::vec intercepts1 = arma::zeros<arma::vec>(K);
  intercepts1(0) = intercepts(0);
  for(int i = 1; i < K; i++){
    intercepts1(i) = intercepts1(i-1) + intercepts(i);
  }
  for(int i = 0; i < K; i++){
    obj_val += linear_coefs(i) * intercepts1(i) + quad_coefs(i) * pow(intercepts1(i),2);
  }
  for(int i = 0; i < (K-1); i++){
    obj_val += std::log(1 - exp(intercepts(i+1))) * logDelta_coefs(i+1);
  }
  obj_val = -obj_val;
  return obj_val;
}


arma::vec optim_ordinal(arma::vec& init_val, arma::vec& linear_coefs,
                        arma::vec& quad_coefs,
                        arma::vec& logDelta_coefs){
  int K = init_val.n_elem;
  // Extract R's optim function
  Rcpp::Environment stats("package:stats");
  Rcpp::Function optim = stats["optim"];
  // Call the optim function from R in C++
  arma::vec uppers = -arma::ones<arma::vec>(K) * 1e-8;
  uppers(0) = 1e4;
  Rcpp::List opt_results = optim(Rcpp::_["par"]    = init_val,
                                 // Make sure this function is not exported!
                                 Rcpp::_["fn"]     = Rcpp::InternalFunction(&obj_fun_ordinal),
                                 Rcpp::_["method"] = "L-BFGS-B",
                                 Rcpp::_["upper"] = uppers,
                                 // Pass in the other parameters as everything
                                 // is scoped environmentally
                                 Rcpp::_["linear_coefs"] = linear_coefs,
                                 Rcpp::_["quad_coefs"] = quad_coefs,
                                 Rcpp::_["logDelta_coefs"] = logDelta_coefs);
  
  // Extract out the estimated parameter values
  arma::vec out = Rcpp::as<arma::vec>(opt_results[0]);
  // Return estimated values
  return out;
}


void ordinalINTERCEPTS(arma::vec& y, arma::vec& yobs, arma::vec& offset, arma::mat& UPXI,
                       int K,  arma::vec& linear_coefs, arma::vec& quad_coefs,
                       arma::vec& logDelta_coefs, const double robust_eps){
  int n = y.n_elem;
  for(int k = 0; k < K; k++){
    arma::uvec idx1 = arma::find((y == k) && (yobs == 1));
    arma::uvec idx2 = arma::find((y == k+1) && (yobs == 1));
    arma::vec lambda = exp(-UPXI.col(k));
    for(int i = 0; i < n; i++){
      if(UPXI(i,k) < 1e-4){
        lambda(i) = 1.0/8.0;
      }else{
        lambda(i) = (1.0 - lambda(i))/(1.0 + lambda(i))/(4 * UPXI(i,k));
      }
    }
    lambda = lambda+robust_eps;
    linear_coefs(k) = (idx2.n_elem * 1.0- idx1.n_elem * 1.0)/2.0;
    linear_coefs(k) -= 2.0 * arma::sum(offset(idx1) % lambda(idx1)) +
      2.0 * arma::sum(offset(idx2) % lambda(idx2));
    quad_coefs(k) = -arma::sum( lambda(idx1)) - arma::sum(lambda(idx2));
    if(k <= K-2){
      logDelta_coefs(k+1) = 1.0 * idx2.n_elem;
    }
  }
}

void update_ordinal_approximation(arma::mat& Y, const arma::mat& Yobs,
                                  const arma::vec nclasses, arma::mat& Yapprox,
                                  arma::mat& meanFactors, arma::mat& U2,
                                  arma::mat& post_tmu, arma::mat& post_tsigma2,
                                  arma::mat& post_tpi, List& intercepts, arma::mat& nu_mat,
                                  const double robust_eps){
  int n = Y.n_rows;
  int p = Y.n_cols;
  int num_factors = meanFactors.n_cols;
  for(int j = 0; j < p; j++){
    //number of classes
    int K = nclasses(j) - 1;
    arma::vec intercept_prev = intercepts(j);
    arma::mat UPXI = arma::zeros<arma::mat>(n, K);
    arma::vec b = (post_tmu.row(j)).t();
    arma::vec b2 = (post_tmu.row(j) % post_tmu.row(j) + post_tsigma2.row(j)).t();
    arma::vec offset = meanFactors * b;
    arma::vec Delta = arma::zeros<arma::vec>(n);
    for(int k = 0; k < num_factors; k++){
      Delta += U2.col(k) * b2(k) - (meanFactors.col(k) % meanFactors.col(k)) * b(k) * b(k);
    }
    for(int k = 0; k < K; k++){
      UPXI.col(k) = (offset + intercept_prev(k));
      UPXI.col(k)  = UPXI.col(k) % UPXI.col(k) + Delta;
      UPXI.col(k) = arma::sqrt(UPXI.col(k));
    }
    //update intercepts
    int m = intercept_prev.n_elem;
    arma::vec init_val = arma::zeros<arma::vec>(m);
    init_val(0) = intercept_prev(0);
    for(int l = 1; l < m; l++){
      init_val(l) = intercept_prev(l)  - intercept_prev(l-1);
      if(init_val(l) >= (-1e-6)){
        init_val(l) = (-1e-6);
      }
    }
    for(int l = 1; l < m; l++){
      intercept_prev(l) = intercept_prev(l-1)  + init_val(l);
    }
    arma::vec y = Y.col(j);
    arma::vec yobs = Yobs.col(j);
    arma::vec linear_coefs = arma::zeros<arma::vec>(K);
    arma::vec quad_coefs = arma::zeros<arma::vec>(K);
    arma::vec logDelta_coefs = arma::zeros<arma::vec>(K);
    ordinalINTERCEPTS(y, yobs, offset, UPXI, K, linear_coefs, quad_coefs, logDelta_coefs, robust_eps);
    arma::vec solve_val = arma::zeros<arma::vec>(m);
    double obj = obj_fun_ordinal(init_val, linear_coefs, quad_coefs, logDelta_coefs);
    solve_val = optim_ordinal(init_val, linear_coefs, quad_coefs,logDelta_coefs);
    intercept_prev(0) = solve_val(0);
    for(int l = 1; l < m; l++){
      intercept_prev(l) = intercept_prev(l-1)  + solve_val(l);
    }
    //update UPXI matrix and Y approx
    UPXI_update_ordinal(j, Y,  Yapprox, offset,  Delta,
                        intercept_prev,  UPXI, nu_mat, robust_eps);
    intercepts(j) = intercept_prev;
  }
}



double binary_search(double val_min, double val_max, double bound, arma::vec diags, arma::vec z, int max_it){
  double cur_val = (val_min + val_max)/2.0;
  arma::vec tmp = 1.0/(diags + cur_val);
  tmp = tmp % tmp;
  double cur_max = val_max;
  double cur_min = val_min;
  double cur_obj = sum(z % tmp);
  double diff = cur_obj - bound;
  int it = 0;
  while((std::abs(diff) > 1e-8) & (it < max_it)){
    if(diff < 0){
      cur_max = cur_val;
      cur_val = (cur_max + cur_min)/2.0;
    }else{
      cur_min = cur_val;
      cur_val = (cur_max + cur_min)/2.0;
    }
    it = it + 1;
    arma::vec tmp1 = 1.0/(diags + cur_val);
    tmp1 = tmp1 % tmp1;
    cur_obj = sum(z % tmp1);
    diff = cur_obj - bound;
  }
  return cur_val;
}


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
 */
double update_projection_sparse(const int num_factors, arma::mat& X,   const arma::mat& Xobs,
                                const arma::vec weights, arma::mat& nu_mat,
                                arma::mat& meanFactors,  arma::mat& U2,
                                arma::mat& post_tmu, arma::mat& post_tsigma2,  arma::mat& post_tpi,
                                arma::mat& tau, arma::mat& log_pi,  arma::mat& log_minus_pi){
  int p = X.n_cols;
  double Delta = 0.0;
  for(int j = 0; j < p; j++){
    arma::uvec obs = arma::find(Xobs.col(j) == 1);
    arma::vec x =  X.col(j);
    x = x(obs);
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
      arma::vec response1 = (response % nu_vec) * weights(j);
      //denominator of eq. (1)
      double tmp = arma::sum(nu_vec % u) * weights(j)+  tau(j,k);
      //numerator of eq. (1)
      double tmp1 = sum(response1 % f);
      //deduct current contribution from B_{jk} to the EBLO using eq (1).
      double delta1 = post_tpi(j,k) * (log_pi(j,k) -log(post_tpi(j,k)));
      double delta2 = (1.0 - post_tpi(j,k)) * (log_minus_pi(j,k) -log(1.0 - post_tpi(j,k)) - 0.5);
      if(post_tpi(j,k) == 0.0){
        delta1 = 0.0;
      }
      if(post_tpi(j,k) == 1.0){
        delta2 = 0.0;
      }
      Delta = Delta - (0.5 * post_tpi(j,k) *(-tmp * bkl2 + 2* tmp1 * post_tmu(j,k) + log(post_tsigma2(j,k) * tau(j,k))) +  delta1+ delta2);
      //update post_tmu, post_tsigma2, post_tpi with formulas (1)-(3)
      post_tsigma2(j,k)  = 1.0/tmp;
      post_tmu(j,k) = tmp1 * post_tsigma2(j,k);
      double tmp2 = 0.5 *(post_tmu(j,k)* post_tmu(j,k)) / post_tsigma2(j,k) +
        0.5 * std::log(tau(j,k) * post_tsigma2(j,k)) +   log_pi(j,k) - log_minus_pi(j,k);
      if(tmp2 < 0){
        post_tpi(j,k) = std::exp(tmp2);
        post_tpi(j,k) = post_tpi(j,k)/(1.0+post_tpi(j,k));
      }else{
        post_tpi(j,k) = std::exp(-tmp2);
        post_tpi(j,k) = 1.0/(1.0+post_tpi(j,k));
      }
      //calculate the new ELBO (partial)
      bkl = post_tmu(j,k) * post_tpi(j,k);
      response = response - f * bkl;
      bkl2 = post_tmu(j,k) * post_tmu(j,k)+post_tsigma2(j,k);
      delta1 = post_tpi(j,k) * (log_pi(j,k) -log(post_tpi(j,k)));
      delta2 = (1.0 - post_tpi(j,k)) * (log_minus_pi(j,k) -log(1.0 - post_tpi(j,k)) - 0.5);
      if(post_tpi(j,k) == 0.0){
        delta1 = 0.0;
      }
      if(post_tpi(j,k) == 1.0){
        delta2 = 0.0;
      }
      Delta = Delta + (0.5 * post_tpi(j,k) *(-tmp * bkl2 + 2* tmp1 * post_tmu(j,k) + log(post_tsigma2(j,k) * tau(j,k)))
                         +  delta1+ delta2);    }
    
  }
  return Delta;
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
*/

double update_projection_constrained(const int num_factors, arma::mat& Y, const arma::mat& Yobs,
                                     arma::mat& nuYmat, arma::mat& meanFactors,  arma::mat& U2,
                                     arma::mat& post_tmu, arma::mat& post_tsigma2,
                                     arma::mat& tau, double lower){
  int p = Y.n_cols;
  double Delta = 0.0;
  for(int j = 0; j < p; j++){
    arma::vec y =  Y.col(j);
    arma::uvec obs = arma::find(Yobs.col(j) == 1);
    arma::vec nu_vec = nuYmat.col(j);
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
          cov(k, k1) = arma::sum(f % f1 % nu_vec);
        }
      }
    }
    for(int k = 0; k < num_factors; k++){
      arma::vec u = U2_sub.col(k);
      cov(k,k) = arma::sum( u % nu_vec)+tau(j,k);
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
    response =meanFactors_sub.t() * (response % nu_vec);
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
 * Calculate E[\|U_k\|_2^2]
 */
arma::mat U2calculation(const int& num_factors, const arma::mat& Z, arma::mat& meanFactors,
                        arma::mat& post_mu, arma::mat& post_sigma2, arma::mat& post_pi){
  int n =meanFactors.n_rows;
  int p = Z.n_cols;
  arma::mat U2 = arma::zeros<arma::mat>(n, num_factors);
  arma::mat diags = arma::zeros<arma::mat>(n, num_factors);
  for(int k = 0; k < num_factors; k++){
    for(int j = 0; j < p; j++){
      double betajk2 = (post_mu(j,k) * post_mu(j,k) +  post_sigma2(j,k)) * post_pi(j,k);
      double betajk = post_mu(j,k)* post_pi(j,k);
      diags.col(k) = diags.col(k) + arma::square(Z.col(j)) * (betajk2 - betajk*betajk);
    }
    U2.col(k)= diags.col(k) +arma::square(meanFactors.col(k));
  }
  return U2;
}


/*
 * update_factor_one: update factors for a given pattern
 * 
 */
double update_factor_one(const int num_factors, arma::mat& Y,  const arma::mat& Yobs,
                         arma::mat& X,  const arma::mat& Xobs,
                         const arma::mat& Z, const arma::vec& weights,
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
      arma::vec yobs = Xobs.col(j);
      arma::uvec obs = arma::find(yobs == 1);
      bkl = post_tmuY(j,k);
      bkl2 = post_tmuY(j,k) * post_tmuY(j,k)+post_tsigma2Y(j,k);
      response(obs) = response(obs)+ (y(obs) %nu_vec(obs)) *bkl;
      s_vec(obs) += bkl2 * nu_vec(obs);
      for(int k1 = 0; k1 < num_factors; k1++){
        arma::vec f = meanFactors.col(k1);
        if(k1 != k){
          double cross = 0.0;
          cross = post_tmuY(j,k1) * post_tmuY(j,k);
          response(obs) = response(obs) - (f(obs) % nu_vec(obs)) * cross;
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
      double tmp =  arma::sum(s_vec % arma::square(Z.col(j)))+tau(j,k);
      //tmp1 = numerator of mu_jk.
      double tmp1 = sum(response1 % Z.col(j));
      double delta1 = post_pi(j,k) * (log_pi(j,k) -log(post_pi(j,k)));
      double delta2 = (1.0 - post_pi(j,k)) * ( log_minus_pi(j,k) -log(1.0 - post_pi(j,k)) - 0.5);
      if(post_pi(j,k) == 0.0){
        delta1 = 0.0;
      }
      if(post_pi(j,k) == 1.0){
        delta2 = 0.0;
      }
      Delta = Delta - (0.5 * post_pi(j,k) *(-tmp * betajk2 + 2* tmp1 * post_mu(j,k) + log(post_sigma2(j,k) * tau(j,k)))
                         +  delta1 + delta2);
      post_sigma2(j,k) = 1.0/tmp;
      post_mu(j,k) = tmp1 *  post_sigma2(j,k);
      double tmp2 = 0.5 * (post_mu(j,k)* post_mu(j,k))/ post_sigma2(j,k) +
        0.5 * std::log(tau(j,k) * post_sigma2(j,k)) + (log_pi(j,k) - log_minus_pi(j,k));
      if(tmp2 < 0){
        post_pi(j,k) = std::exp(tmp2);
        post_pi(j,k) = post_pi(j,k)/(1.0+post_pi(j,k));
      }else{
        post_pi(j,k) = std::exp(-tmp2);
        post_pi(j,k) = 1.0/(1.0+post_pi(j,k));
      }
      betajk2 = post_mu(j,k) * post_mu(j,k)+post_sigma2(j,k);
      betajk = post_mu(j,k) * post_pi(j,k);
      delta1 = post_pi(j,k) * (log_pi(j,k) -log(post_pi(j,k)));
      delta2 = (1.0 - post_pi(j,k)) * ( log_minus_pi(j,k) -log(1.0 - post_pi(j,k)) - 0.5);
      if(post_pi(j,k) == 0.0){
        delta1 = 0.0;
      }
      if(post_pi(j,k) == 1.0){
        delta2 = 0.0;
      }
      Delta = Delta + (0.5 * post_pi(j,k) *(-tmp * betajk2 + 2* tmp1 * post_mu(j,k) + log(post_sigma2(j,k) * tau(j,k)))
                         +  delta1 + delta2);
      meanFactors.col(k) += Z.col(j) * betajk;
    }
  }
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
 * 
 */
double update_factor(const int num_factors, arma::mat& Y,  const arma::mat& Yobs,
                     arma::mat& X,  const arma::mat& Xobs,
                     const arma::mat& Z,  const arma::vec& weights,
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
                                      Zsub, weights, post_mu_sub, post_sigma2_sub,post_pi_sub ,
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
 * Update the prior distribution of coefficient variance for inactive features.
 */
double tau_update(const int& num_factors, const arma::vec& weights,
                  const List& pattern_features, const List& functional_path,
                  const double a0, const double b0,arma::cube& post_mu,
                  arma::cube& post_sigma2, arma::cube& post_pi,
                  arma::mat& post_tmu, arma::mat& post_tsigma2,arma::mat& post_tpi,
                  arma::mat& tauZ, arma::mat& post_a0, arma::mat& post_b0, const double L){
  //posterior distribution for tau
  int num_path = functional_path.size();
  int num_pattern = pattern_features.size();
  double Delta = 0.0;
  for(int l = 0; l < num_path; l++){
    arma::uvec path_index = functional_path(l);
    int p0 = path_index.n_elem;
    for(int k = 0; k < num_factors; k++){
      double post_a = a0;
      double post_b = b0;
      //add up the effects from factors
      for(int q = 0; q < num_pattern; q++){
        arma::uvec jj = pattern_features(q);
        arma::uvec path_q = intersect(jj, path_index);
        int p1 = path_q.n_elem;
        for(int j0 = 0; j0 < p1; j0++){
          int j = path_q(j0);
          double tmp = (post_mu(j, k, q) * post_mu(j, k, q) + post_sigma2(j,k,q)) * post_pi(j,k,q) * .5;
          post_a = post_a + post_pi(j,k,q) * .5;
          post_b = post_b + tmp;
        }
      }
      //add up the effects from weights
      for(int j0 = 0; j0 < p0; j0++){
        int j = path_index(j0);
        if(weights(j) < 1){
          post_a = post_a +   post_tpi(j,k) * .5* weights(j);
          post_b = post_b +   ((post_tmu(j, k) * post_tmu(j, k) + post_tsigma2(j,k)) * post_tpi(j,k)) * .5* weights(j);
        }else{
          post_a = post_a +   post_tpi(j,k) * .5;
          post_b = post_b +   ((post_tmu(j, k) * post_tmu(j, k) + post_tsigma2(j,k)) * post_tpi(j,k)) * .5;
        }
      }
      Delta = Delta - ((post_a-post_a0(l,k)) * (boost::math::digamma(post_a0(l,k))-log(post_b0(l,k)))
                         - (post_b-post_b0(l,k)) * post_a0(l,k)/post_b0(l,k) +
                           std::lgamma(post_a0(l,k)) - post_a0(l,k) * log(post_b0(l,k)));
      double L_est = post_a/post_b;
      double post_a_use = post_a;
      double post_b_use = 0.0;
      if(L_est > L){
        post_b_use = post_a_use/L;
      }else{
        post_b_use = post_b;
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
  return Delta;
}


/*
 * Update the prior distribution of sparsity levels.
 */
double pi_update(const int& num_factors, const arma::vec& weights,
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
    int p0 = path_index.n_elem;
    for(int k = 0; k < num_factors; k++){
      double weight_sum = 0.0;
      double total_sum = 0.0;
      //add up the effects from factors
      for(int q = 0; q < num_pattern; q++){
        arma::uvec jj = pattern_features(q);
        arma::uvec path_q = intersect(jj, path_index);
        int p1 = path_q.n_elem;
        for(int j0 = 0; j0 < p1; j0++){
          int j = path_q(j0);
          total_sum = total_sum + 1.0;
          weight_sum = weight_sum + post_pi(j,k,q);
        }
      }
      //add up the effects from weights
      for(int j0 = 0; j0 < p0; j0++){
        int j = path_index(j0);
        if(weights(j) < 1.0){
          weight_sum = weight_sum +post_tpi(j,k) * weights(j);
          total_sum = total_sum + weights(j);
        }else{
          weight_sum = weight_sum +post_tpi(j,k);
          total_sum = total_sum + 1.0;
        }
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
      if(post_a1(l,k)>alpha0 * (post_a1(l,k)+post_b1(l,k))){
        post_a1(l,k) = alpha0 * (total_sum+b1+a1);
        post_b1(l,k) = (1.0-alpha0) * (total_sum+b1+a1);
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
 */
double nu_update(arma::mat& Y, const arma::mat& Yobs, const arma::vec& weights,
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
      LambdaU(k,k) = arma::sum(u2(obs))-arma::sum(u20(obs));
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
      correction += LambdaB(k,k) * arma::sum(u20(obs));
      correction += tmp2 * LambdaU(k,k);
      //correction = correction + arma::sum(u2(obs)) * tmp2 - tmp1*tmp1*arma::sum(u20(obs));
    }
    if(!sparse){
      //correction = 0.0;
    }
    double tmp =  arma::sum(arma::square(r(obs))) + correction;
    double post_a = a2 + weights(j) * n/2.0;
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

