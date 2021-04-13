#ifndef PKG_SPEAR01_H
#define PKG_SPEAR01_H
#include <RcppArmadillo.h>
#include <boost/math/special_functions/digamma.hpp>
#include <cmath>
// [[Rcpp::depends(BH)]]
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

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
                                arma::mat& tau, arma::mat& log_pi,  arma::mat& log_minus_pi);

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
                                     arma::mat& tau, double lower);

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
                     arma::mat& updatingOrders);

double update_factor_one(const int num_factors, arma::mat& Y,  const arma::mat& Yobs,
                         arma::mat& X,  const arma::mat& Xobs,
                         const arma::mat& Z, const arma::vec& weights,
                         arma::mat& post_mu, arma::mat& post_sigma2, arma::mat& post_pi,
                         arma::mat& post_tmuX, arma::mat& post_tsigma2X,
                         arma::mat& post_tpiX, arma::mat& post_tmuY, arma::mat& post_tsigma2Y,
                         arma::mat& tau,  arma::mat& log_pi, arma::mat& log_minus_pi,
                         arma::mat& meanFactors, arma::mat& nuXmat,
                         arma::mat& nuYmat, arma::mat& updatingOrders);

double tau_update(const int& num_factors, const arma::vec& weights,
                  const List& pattern_features, const List& functional_path,
                  const double a0, const double b0,arma::cube& post_mu,
                  arma::cube& post_sigma2, arma::cube& post_pi,
                  arma::mat& post_tmu, arma::mat& post_tsigma2,arma::mat& post_tpi,
                  arma::mat& tauZ, arma::mat& post_a0, arma::mat& post_b0, const double L);

double pi_update(const int& num_factors, const arma::vec& weights,
                 const List& pattern_features, const List& functional_path,
                 const double a1, const double b1,
                 arma::cube& post_mu, arma::cube& post_sigma2, arma::cube& post_pi,
                 arma::mat& post_tmu, arma::mat& post_tsigma2,arma::mat& post_tpi,
                 arma::mat& log_pi, arma::mat& log_minus_pi,arma::mat& post_a1, arma::mat& post_b1, const double alpha0);

double nu_update(arma::mat& Y, const arma::mat& Yobs, const arma::vec& weights,
                 arma::mat& meanFactors, arma::mat& U2,arma::mat& nu_mat,
                 const double a2, const double b2,
                 arma::mat& post_tmu, arma::mat& post_tsigma2,arma::mat& post_tpi,
                 arma::vec& post_a2, arma::vec& post_b2, bool sparse);


double binary_search(double val_min, double val_max, double bound, arma::vec diags, arma::vec z, int max_it);


arma::mat U2calculation(const int& num_factors, const arma::mat& Z, arma::mat& meanFactors,
                        arma::mat& post_mu, arma::mat& post_sigma2, arma::mat& post_pi);

void UPXI_update_binomial(int idx, arma::mat& S, arma::mat& Y, arma::mat& Yapprox,
                          arma::vec& offset, arma::vec& Delta,
                          arma::vec& intercepts, arma::mat& UPXI, arma::mat& nu_mat,
                          const double robust_eps);

void update_binomial_approximation(arma::mat& Y, const arma::mat& Yobs,
                                   const arma::vec nclasses,  arma::mat& Yapprox,
                                   arma::mat& meanFactors, arma::mat& U2,
                                   arma::mat& post_tmu, arma::mat& post_tsigma2,
                                   arma::mat& post_tpi, List& intercepts, arma::mat& nu_mat,
                                   const double robust_eps);

void UPXI_update_ordinal(int idx,  arma::mat& Y,  arma::mat& Yapprox,
                         arma::vec& offset, arma::vec& Delta,
                         arma::vec& intercepts, arma::mat& UPXI, arma::mat& nu_mat,
                         const double robust_eps);

double obj_fun_ordinal(arma::vec& intercepts,
                       arma::vec& linear_coefs, arma::vec& quad_coefs, arma::vec& logDelta_coefs);

arma::vec optim_ordinal(arma::vec& init_val, arma::vec& linear_coefs,
                        arma::vec& quad_coefs,
                        arma::vec& logDelta_coefs);

void ordinalINTERCEPTS(arma::vec& y, arma::vec& yobs, arma::vec& offset, arma::mat& UPXI,
                       int K,  arma::vec& linear_coefs, arma::vec& quad_coefs,
                       arma::vec& logDelta_coefs, const double robust_eps);


void update_ordinal_approximation(arma::mat& Y, const arma::mat& Yobs,
                                  const arma::vec nclasses, arma::mat& Yapprox,
                                  arma::mat& meanFactors, arma::mat& U2,
                                  arma::mat& post_tmu, arma::mat& post_tsigma2,
                                  arma::mat& post_tpi, List& intercepts, arma::mat& nu_mat,
                                  const double robust_eps);

double binary_search(double val_min, double val_max, double bound, arma::vec diags, arma::vec z, int max_it);

#endif
