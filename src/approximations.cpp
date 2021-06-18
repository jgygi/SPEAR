#include <RcppArmadillo.h>
#include <boost/math/special_functions/digamma.hpp>
#include <cmath>
// [[Rcpp::depends(BH)]]
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

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
  arma::vec uppers = -arma::ones<arma::vec>(K) * 1e-6;
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


/*
 * Logistic multinomial
 */
void UPXI_update_multinomial(arma::mat& Y, arma::mat& Yapprox,
                             arma::mat& offset, arma::mat& Delta,
                             arma::vec& intercepts, arma::mat& UPXI,  arma::vec& UPXIjoint, arma::mat& nu_mat,
                             const double robust_eps){
  //Boucher bound
  int n = UPXI.n_rows;
  int p = UPXI.n_cols;
  //iterate three times
  for(int j = 0; j < p; j++){
    UPXI.col(j) = (offset.col(j) + intercepts(j)-UPXIjoint);
    UPXI.col(j)  = UPXI.col(j) % UPXI.col(j) + Delta.col(j);
    UPXI.col(j) = arma::sqrt(UPXI.col(j));
  }
  UPXIjoint = arma::ones<arma::vec>(n) * 1/2.0 * (p/2.0 - 1);
  arma::vec tmp0 = arma::zeros<arma::vec>(n);
  for(int j = 0; j < p; j++){
    arma::vec lambdas = exp(-UPXI.col(j));
    lambdas = (1.0 - lambdas)/(1.0+lambdas)/(4 * UPXI.col(j))+robust_eps/p;
    UPXIjoint += lambdas % (offset.col(j)+intercepts(j));
    tmp0 += lambdas;
  }
  UPXIjoint  = UPXIjoint /tmp0;
}

void Yapprox_update_multinomial(arma::mat& Y, arma::mat& Yapprox,
                                arma::mat& offset, arma::mat& Delta,
                                arma::vec& intercepts, arma::mat& UPXI,  arma::vec& UPXIjoint, arma::mat& nu_mat,
                                const double robust_eps){
  //Boucher bound
  int n = UPXI.n_rows;
  int p = UPXI.n_cols;
  for(int j = 0; j < p; j++){
    arma::vec lambdas = exp(-UPXI.col(j));
    lambdas = (1.0 - lambdas)/(1.0+lambdas)/(4 * UPXI.col(j));
    for(int i = 0; i < n; i++){
      if(UPXI(i,j) < 1e-4){
        lambdas(i) = 1.0/8.0;
      }
    }
    lambdas = lambdas + robust_eps/p;
    nu_mat.col(j) = 2.0 * lambdas;
    Yapprox.col(j) = (Y.col(j) - 0.5  - 2* lambdas % UPXIjoint)/nu_mat.col(j)  -  intercepts(j);
  }
}

void update_multinomial_approximation(arma::mat& Y, const arma::mat& Yobs,
                                      arma::mat& UPXI,  arma::vec& UPXIjoint,
                                      const arma::vec nclasses,  arma::mat& Yapprox,
                                      arma::mat& meanFactors, arma::mat& U2,
                                      arma::mat& post_tmu, arma::mat& post_tsigma2,
                                      arma::mat& post_tpi, List& intercepts, arma::mat& nu_mat,
                                      const double robust_eps){
  int p = Y.n_cols;
  int n = Y.n_rows;
  int num_factors = meanFactors.n_cols;
  arma::mat offset = arma::zeros<arma::mat>(n,p);
  arma::mat Delta = arma::zeros<arma::mat>(n,p);
  for(int j = 0; j < p; j++){
    //number of classes
    arma::vec intercept_prev = intercepts(j);
    arma::vec b = (post_tmu.row(j) % post_tpi.row(j)).t();
    arma::vec b2 = ((post_tmu.row(j) % post_tmu.row(j) + post_tsigma2.row(j)) % post_tpi.row(j)).t();
    offset.col(j) = meanFactors * b;
    for(int k = 0; k < num_factors; k++){
      Delta.col(j) += U2.col(k) * b2(k) - (meanFactors.col(k) % meanFactors.col(k)) * b(k) * b(k);
    }
  }
  arma::vec tmp_intercepts = arma::zeros<arma::vec>(p);
  for(int j = 0; j < p; j++){
    tmp_intercepts(j) = intercepts(j);
  }
  UPXI_update_multinomial(Y,  Yapprox, offset, Delta, tmp_intercepts, UPXI, UPXIjoint, nu_mat,robust_eps);
  //update intercepts USING the non-missing values!
  for(int j = 0; j < p; j++){
    arma::vec lambdas = exp(-UPXI.col(j));
    lambdas = (1.0 - lambdas)/(1.0+lambdas)/(4 * UPXI.col(j))+robust_eps/p;
    arma::vec tmp1 = (Y.col(j) - 0.5 - 2 * lambdas % UPXIjoint)/(2 * lambdas);
    tmp1 = tmp1 - offset.col(j);
    tmp_intercepts(j)= arma::sum(tmp1 % lambdas)/arma::sum(lambdas);
    intercepts(j) =tmp_intercepts(j); 
  }
  //update UPXI matrix and Y approx
  Yapprox_update_multinomial(Y, Yapprox, offset, Delta, tmp_intercepts, UPXI,  UPXIjoint, nu_mat,robust_eps);
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
  while((std::abs(diff) > 1e-12) & (it < max_it)){
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
