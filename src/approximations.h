#ifndef PKG_SPEAR02_H
#define PKG_SPEAR02_H
#include <RcppArmadillo.h>
#include <boost/math/special_functions/digamma.hpp>
#include <cmath>
// [[Rcpp::depends(BH)]]
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

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

void update_multinomial_approximation(arma::mat& Y, const arma::mat& Yobs,
                                      arma::mat& UPXI,  arma::vec& UPXIjoint,
                                      const arma::vec nclasses,  arma::mat& Yapprox,
                                      arma::mat& meanFactors, arma::mat& U2,
                                      arma::mat& post_tmu, arma::mat& post_tsigma2,
                                      arma::mat& post_tpi, List& intercepts, arma::mat& nu_mat,
                                      const double robust_eps);


double binary_search(double val_min, double val_max, double bound, arma::vec diags, arma::vec z, int max_it);


arma::mat U2calculation(const int& num_factors, const arma::mat& Z, arma::mat& meanFactors,
                        arma::mat& post_mu, arma::mat& post_sigma2, arma::mat& post_pi);

#endif
