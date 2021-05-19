#include "subroutines.h"
inline int randWrapper(const int n){return floor(unif_rand()*n);}

void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

double pos(double x){
  if(x <= 0.0){
    x = 0.0;
  }
  return x;
}


/* spear_ is the function for external use
  family: datatype of y, 0 = Gaussian, 1 = binomial, 2 = ordinal, 3 = multinomial
  Y: N * p0 response matrix.
  X: N* p feature omics matrix for it.
  Yobs: N*p0 indicator matrix of wether Y is observed, Yobs[i,j] = 1 if Y[i,j] is observed.
  Xobs: N*p indicator matrix of wether X is observed, Xobs[i,j] = 1 if X[i,j] is observed.
  Z: N * p1 (p1>=p) feature matrix for constructing factors U = Zbeta. Z does not contain missing data, one common choice is to let Z be the imputed version of X.
  nclasses: Lenth p0 vector indicating number of unique classes in Y for ordinal regression.
  functional_path: lists of features that should be grouped together.
  pattern_samples: lists of samples having a common beta for each sample group.
  pattern_features: lists of features used in beta for each sample group.
  weights: weights for X.
  thres_eblo, thres_count: stop if increase of eblo is below thre_eblo for consequtive thres_count.
  thres_factor: drop a factor if its variance is below thres_factor.
  (a0, b0), (a1, b1), (a2, b2): hyper priors for sparsity, effect size and noise variance.
  lower: ??
  tauY, tauZ: variance of B for Y, X/Z (?)
  log_pi, log_pi_minus: expected log pi and log(1-pi) of sparsity level for beta and B.
  (post_a0, post_b0): posteriors Gamma distribution paramters for sparisy levels
  (post_a1, post_b1): posteriors Beta distribution paramters for beta
  (post_ax, post_bx): posteriors Beta distribution paramters for B in X
  (post_ay, post_by): posteriors Beta distribution paramters for B in Y
  robust_eps: small add-on quantities for nonGaussian data to avoid numerical error in the completely separable case.
  alpha0: sparisty upper bound on the sparsity level per pathway.
  # outputs directly used for downstream analysis
  interceptsX, interceptsY: intercepts
  post_mu, post_sigma2, post_pi: posterior distribution of beta when non-zero and probability of being non-zero.
  post_tmuX, post_tsigma2X, post_tpiX: posterior distribution of B in X when non-zero and probability of being non-zero.
  post_tmuY, post_tsigma2Y, post_tpiY: posterior distribution of B in Y when non-zero and probability of being non-zero.
 post_tpiY is not active since there is no selection on B for Y.
  
 */

//' @export
// [[Rcpp::export]]
arma::mat spear_(const int family, arma::mat& Y,  arma::mat& X,
            const arma::mat& Yobs, const arma::mat& Xobs,
            const arma::mat& Z, const arma::vec nclasses, const List& functional_path,
            const List& pattern_samples, const List& pattern_features,
            const arma::vec weights, const int num_factors, 
            const int warm_up, const int max_iter,
            const double thres_elbo, const int thres_count, const double thres_factor, 
            const double a0, const double b0, const double a1, const double b1,
            const double a2, const double b2, const double lower, const int print_out,
            arma::vec interceptsX, List& interceptsY, 
            arma::cube& post_mu, arma::cube& post_sigma2, arma::cube& post_pi, 
            arma::mat& post_tmuX, arma::mat& post_tsigma2X, arma::mat& post_tpiX, 
            arma::mat& post_tpiX_marginal, arma::mat& post_tmuY, arma::mat& post_tsigma2Y, arma::mat& post_tpiY,
            arma::mat& tauY, arma::mat& tauZ,
            arma::mat& log_pi,arma::mat& log_minus_pi, 
            arma::mat& nuXmat, arma::mat& nuYmat,
            arma::mat& post_a0, arma::mat& post_b0,
            arma::mat& post_a1, arma::mat& post_b1,
            arma::vec& post_a2x, arma::vec& post_b2x, 
            arma::vec& post_a2y, arma::vec& post_b2y,
            arma::mat& meanFactors, const int seed0,
            const double robust_eps, const double alpha0, const double L){
    int num_pattern = pattern_samples.size();
    int py = Y.n_cols;
    int px = X.n_cols;
    int pz = Z.n_cols;
    int n = Y.n_rows;
    int consecutives_count = 0;
    arma::mat UPXI = arma::zeros<arma::mat>(n, py);
    arma::vec UPXIjoint = arma::zeros<arma::vec>(n);
    /*
    Yapprox and Xapprox can be viewed as residuals that we want to further fit to.
    It will be updated based on current model fits. Update can occur due to
    updated model fit and as well as updated approximation for non-Gaussian Y.
    */
    arma::mat Yapprox = Y * 1.0;
    arma::mat Xapprox = X * 1.0;
    /*
    Expected values of the square factor value.
   */
    arma::mat U2 = arma::zeros<arma::mat>(n, num_factors);
    /*
     delta* is the increase in EBLO for each subroutines and Delta is the total increase.
    */
    double delta11, delta12, delta2, delta3, delta4, delta51, delta52 = 0.0;
    double Delta = 0.0;
    /*
     For easy control of reproducibility, we will contain all randomness at this chunck
     of code. This including:
     1) generate all updating orders for later iterations.
     */
    set_seed(seed0);
    arma::cube updatingOrdersAll = arma::zeros<arma::cube>(pz, num_factors, max_iter);
    for(int it = 0; it < max_iter; it++){
        for(int k = 0; k < num_factors; k++){
            IntegerVector shuffles_feature = seq(0,pz-1);
            std::random_shuffle(shuffles_feature.begin(),shuffles_feature.end(),randWrapper);
            for(int j = 0; j < pz; j++){
                updatingOrdersAll(j,k,it) = shuffles_feature(j);
            }
        }
    }
    /*
     We update model parameters in the following orders
     1) Initialization
     1) B from factor to resonse
     2) beta -> factor: Update for pathway separately given priors and B.
     3) B from factor to X
     */
    //Initialization, step 0: calculate the expected factor and expected square factor given posterior distribution of beta. This is done for each pattern separately.
    for(int k = 0; k < num_pattern; k++){
        arma::uvec ii = pattern_samples(k);
        arma::uvec jj = pattern_features(k);
        ii = ii - 1;
        jj = jj - 1;
        const arma::mat Zsub = Z.submat(ii,jj);
        arma::mat post_mu_sub = post_mu.slice(k);
        arma::mat post_sigma2_sub = post_sigma2.slice(k);
        arma::mat post_pi_sub = post_pi.slice(k);
        post_mu_sub = post_mu_sub.rows(jj);
        post_sigma2_sub = post_sigma2_sub.rows(jj);
        post_pi_sub = post_pi_sub.rows(jj);
        //calculate the expected factor values
        arma::mat meanFactors_sub = Zsub * (post_mu_sub % post_pi_sub);
        //calculate the expected squared factor values
        arma::mat U2_sub =  U2calculation(num_factors, Zsub, meanFactors_sub,
                                          post_mu_sub, post_sigma2_sub, post_pi_sub);
        U2.rows(ii) = U2_sub;
        meanFactors.rows(ii) = meanFactors_sub;
    }
    // Initialization step1:
    //remove intercepts
    if(family == 0){
      for(int j = 0; j < py; j++){
        arma::uvec obs = arma::find(Yobs.col(j) == 1);
        arma::vec y = Y.col(j) * 1.0;
        arma::vec b = (post_tmuY.row(j)).t();
        arma::vec yhat = meanFactors * b;
        interceptsY(j)  = (arma::mean(y(obs) - yhat(obs))) *arma::ones<arma::vec>(1);
        arma::vec intercept0 = interceptsY(j);
        Yapprox.col(j) = y - intercept0(0);
      }
    }
    for(int j = 0; j < px; j++){
      arma::uvec obs = arma::find(Xobs.col(j) == 1);
      arma::vec x = X.col(j) * 1.0;
      arma::vec b = (post_tmuX.row(j) % post_tpiX.row(j)).t();
      arma::vec xhat = meanFactors * b;
      interceptsX(j) = arma::mean(x(obs) - xhat(obs));
      Xapprox.col(j) = x - interceptsX(j);
    }
    //Update the residuals given the current model.
    if(family == 1){
        update_binomial_approximation(Y, Yobs, nclasses, Yapprox, meanFactors,
                                      U2,  post_tmuY, post_tsigma2Y,  post_tpiY, interceptsY, nuYmat, robust_eps);
    }else if(family == 2){
        update_ordinal_approximation(Y, Yobs, nclasses, Yapprox, meanFactors,
                                     U2, post_tmuY, post_tsigma2Y,  post_tpiY,  interceptsY, nuYmat, robust_eps);
    }else if(family == 3){
        update_multinomial_approximation(Y, Yobs, UPXI, UPXIjoint, nclasses,Yapprox, meanFactors, U2,
                                       post_tmuY, post_tsigma2Y, post_tpiY, interceptsY, nuYmat, robust_eps);
    }
   //Initiliazation done. Start iteration
    int it = 0;
    while(((it < max_iter) & (consecutives_count < thres_count)) | (it < warm_up)){
        // Update posteriors of beta for X and Y, Update approximated response values and variance for non-Gaussian data.
        arma::mat updatingOrders = updatingOrdersAll.slice(it);
        delta11 = update_projection_sparse(num_factors, Xapprox, Xobs, weights, nuXmat, meanFactors, U2,
                                           post_tmuX, post_tsigma2X, post_tpiX, tauZ,
                                           log_pi, log_minus_pi);
        delta12 = update_projection_constrained(num_factors, Yapprox, Yobs, nuYmat, meanFactors, U2,
                                                post_tmuY,  post_tsigma2Y,
                                                tauY,  lower);
        // Update for B in constructing factors
        delta2 = update_factor(num_factors, Yapprox,  Yobs, Xapprox,  Xobs, Z,  weights,
                               pattern_samples, pattern_features,
                               post_mu, post_sigma2, post_pi,
                               post_tmuX, post_tsigma2X, post_tpiX,
                               post_tmuY, post_tsigma2Y,
                               tauZ,  log_pi, log_minus_pi,
                               nuXmat, nuYmat, meanFactors,  U2, updatingOrders);
        //Update the residuals given the current model.
        for(int j = 0; j < px; j++){
          arma::uvec obs = arma::find(Xobs.col(j) == 1);
          arma::vec x = X.col(j) * 1.0;
          arma::vec b = (post_tmuX.row(j)  % post_tpiX.row(j)).t();
          arma::vec xhat = meanFactors * b;
          interceptsX(j) = arma::mean(x(obs) - xhat(obs));
          Xapprox.col(j) =  x - interceptsX(j);
        }
        if(family == 0){
          for(int j = 0; j < py; j++){
            arma::uvec obs = arma::find(Yobs.col(j) == 1);
            arma::vec y = Y.col(j);
            arma::vec b = (post_tmuY.row(j)).t();
            arma::vec yhat = meanFactors * b;
            interceptsY(j) = (arma::mean(y(obs) - yhat(obs))) *arma::ones<arma::vec>(1);
            arma::vec intercept0 = interceptsY(j);
            Yapprox.col(j) = y - intercept0(0);
          }
        }else if(family == 1){
          update_binomial_approximation(Y, Yobs, nclasses, Yapprox, meanFactors,
                                        U2,  post_tmuY, post_tsigma2Y,  post_tpiY, interceptsY, nuYmat, robust_eps);
        }else if(family == 2){
          update_ordinal_approximation(Y, Yobs, nclasses, Yapprox, meanFactors,
                                       U2, post_tmuY, post_tsigma2Y,  post_tpiY,  interceptsY, nuYmat, robust_eps);
        }else if(family == 3){
          update_multinomial_approximation(Y, Yobs, UPXI, UPXIjoint, nclasses,Yapprox, meanFactors, U2,
                                           post_tmuY, post_tsigma2Y, post_tpiY, interceptsY, nuYmat, robust_eps);
        }
        //update prior tau
        delta3 = tau_update(num_factors, weights, pattern_features , functional_path,
                            a0, b0 , post_mu, post_sigma2, post_pi,
                            post_tmuX, post_tsigma2X, post_tpiX,
                            tauZ, post_a0, post_b0,L);
        //update prior pi
        delta4 = pi_update(num_factors, weights, pattern_features ,functional_path,
                           a1, b1, post_mu, post_sigma2, post_pi,
                           post_tmuX, post_tsigma2X, post_tpiX,
                           log_pi, log_minus_pi, post_a1, post_b1, alpha0);
        //update prior nu
        delta51 = nu_update(Xapprox, Xobs, weights, meanFactors, U2, nuXmat, a2, b2,
                            post_tmuX, post_tsigma2X, post_tpiX,post_a2x, post_b2x, true);
        if(family == 0){
            delta52 = nu_update(Yapprox, Yobs, weights, meanFactors, U2, nuYmat, a2, b2,
                           post_tmuY, post_tsigma2Y, post_tpiY,post_a2y, post_b2y, false);
        }
        delta11 = pos(delta11);
        delta12 = pos(delta12);
        delta2 = pos(delta2);
        delta3 = pos(delta3);
        delta4 = pos(delta4);
        delta51 = pos(delta51);
        delta52 = pos(delta52);
        Delta = delta11 + delta2 + delta3 + delta4 + delta51 + delta52;
        if((Delta  < thres_elbo) & (it >= warm_up)){
            consecutives_count = consecutives_count + 1;
        }else{
            consecutives_count = 0;
        }
        if(it <warm_up){
            for(int k = 0; k < num_factors; k++){
                for(int j = 0; j < px; j++){
                    log_pi(j,k) = sqrt(it/warm_up) * log_pi(j,k) + sqrt((warm_up - it)/warm_up) * std::log(.5);
                    log_minus_pi(j,k) = sqrt(it/warm_up) * log_minus_pi(j,k) + sqrt((warm_up - it)/warm_up) * std::log(.5);
                }
            }
        }
        if((it > warm_up) & ((it - warm_up)% print_out == 0)){
            Rcout << "###############iteration" << it -warm_up <<
            "############### EBLO increase" << Delta <<  "###############"<< "\n";
        }
        it += 1;
    }
    it = 0;
    while(it < 2){
        arma::vec ones = arma::ones<arma::vec>(px);
        delta11 = update_projection_sparse(num_factors, Xapprox, Xobs, ones, nuXmat, meanFactors, U2,
                                           post_tmuX, post_tsigma2X, post_tpiX, tauZ,
                                           log_pi, log_minus_pi);
        delta51 =  nu_update(Xapprox, Xobs, ones, meanFactors, U2, nuXmat, a2, b2,
                             post_tmuX, post_tsigma2X, post_tpiX,post_a2x, post_b2x, true);
        it = it + 1;
    }
    //calculate the marginal significance
    for(int k = 0; k < num_factors; k++){
      arma::vec meanfac0 = meanFactors.col(k);
      arma::vec u20 = U2.col(k);
      for(int j = 0; j < px; j++){
        arma::vec nu_vec = nuXmat.col(j);
        arma::vec x = Xapprox.col(j);
        arma::uvec obs = arma::find(Xobs.col(j) == 1);
        x = x(obs);
        nu_vec = nu_vec(obs);
        arma::vec meanfac = meanfac0(obs);
        arma::vec u2 = u20(obs);
        double taux = tauZ(j,k);
        double log_pix = log_pi(j,k);
        double log_minus_pix = log_minus_pi(j,k);
        post_tpiX_marginal(j,k) = update_marginal_prob(x,  nu_vec, meanfac, u2, taux, log_pix, log_minus_pix);
      }
    }
    return Yapprox;
}



