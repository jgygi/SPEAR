// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// spear_
arma::mat spear_(const int family, arma::mat& Y, arma::mat& X, const arma::mat& Yobs, const arma::mat& Xobs, const arma::mat& Z, const arma::vec nclasses, const List& functional_path, const List& pattern_samples, const List& pattern_features, const arma::vec weights, const arma::vec weights0, const arma::vec weights_case, const int num_factors, const int warm_up, const int max_iter, const double thres_elbo, const int thres_count, const double thres_factor, const double a0, const double b0, const double a1, const double b1, const double a2, const double b2, const double lower, const int print_out, arma::vec interceptsX, List& interceptsY, arma::cube& post_mu, arma::cube& post_sigma2, arma::cube& post_pi, arma::mat& post_tmuX, arma::mat& post_tsigma2X, arma::mat& post_tpiX, arma::mat& post_tpiX_marginal, arma::mat& post_tmuY, arma::mat& post_tsigma2Y, arma::mat& post_tpiY, arma::mat& tauY, arma::mat& tauZ, arma::mat& log_pi, arma::mat& log_minus_pi, arma::mat& nuXmat, arma::mat& nuYmat, arma::mat& post_a0, arma::mat& post_b0, arma::mat& post_a1, arma::mat& post_b1, arma::vec& post_a2x, arma::vec& post_b2x, arma::vec& post_a2y, arma::vec& post_b2y, arma::mat& meanFactors, const int seed0, const double robust_eps, const double alpha0, const double L, const double L2);
RcppExport SEXP _SPEAR_spear_(SEXP familySEXP, SEXP YSEXP, SEXP XSEXP, SEXP YobsSEXP, SEXP XobsSEXP, SEXP ZSEXP, SEXP nclassesSEXP, SEXP functional_pathSEXP, SEXP pattern_samplesSEXP, SEXP pattern_featuresSEXP, SEXP weightsSEXP, SEXP weights0SEXP, SEXP weights_caseSEXP, SEXP num_factorsSEXP, SEXP warm_upSEXP, SEXP max_iterSEXP, SEXP thres_elboSEXP, SEXP thres_countSEXP, SEXP thres_factorSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP a1SEXP, SEXP b1SEXP, SEXP a2SEXP, SEXP b2SEXP, SEXP lowerSEXP, SEXP print_outSEXP, SEXP interceptsXSEXP, SEXP interceptsYSEXP, SEXP post_muSEXP, SEXP post_sigma2SEXP, SEXP post_piSEXP, SEXP post_tmuXSEXP, SEXP post_tsigma2XSEXP, SEXP post_tpiXSEXP, SEXP post_tpiX_marginalSEXP, SEXP post_tmuYSEXP, SEXP post_tsigma2YSEXP, SEXP post_tpiYSEXP, SEXP tauYSEXP, SEXP tauZSEXP, SEXP log_piSEXP, SEXP log_minus_piSEXP, SEXP nuXmatSEXP, SEXP nuYmatSEXP, SEXP post_a0SEXP, SEXP post_b0SEXP, SEXP post_a1SEXP, SEXP post_b1SEXP, SEXP post_a2xSEXP, SEXP post_b2xSEXP, SEXP post_a2ySEXP, SEXP post_b2ySEXP, SEXP meanFactorsSEXP, SEXP seed0SEXP, SEXP robust_epsSEXP, SEXP alpha0SEXP, SEXP LSEXP, SEXP L2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type family(familySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Yobs(YobsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Xobs(XobsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type nclasses(nclassesSEXP);
    Rcpp::traits::input_parameter< const List& >::type functional_path(functional_pathSEXP);
    Rcpp::traits::input_parameter< const List& >::type pattern_samples(pattern_samplesSEXP);
    Rcpp::traits::input_parameter< const List& >::type pattern_features(pattern_featuresSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type weights0(weights0SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type weights_case(weights_caseSEXP);
    Rcpp::traits::input_parameter< const int >::type num_factors(num_factorsSEXP);
    Rcpp::traits::input_parameter< const int >::type warm_up(warm_upSEXP);
    Rcpp::traits::input_parameter< const int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const double >::type thres_elbo(thres_elboSEXP);
    Rcpp::traits::input_parameter< const int >::type thres_count(thres_countSEXP);
    Rcpp::traits::input_parameter< const double >::type thres_factor(thres_factorSEXP);
    Rcpp::traits::input_parameter< const double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< const double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< const double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< const double >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< const double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< const double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const int >::type print_out(print_outSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type interceptsX(interceptsXSEXP);
    Rcpp::traits::input_parameter< List& >::type interceptsY(interceptsYSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type post_mu(post_muSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type post_sigma2(post_sigma2SEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type post_pi(post_piSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type post_tmuX(post_tmuXSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type post_tsigma2X(post_tsigma2XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type post_tpiX(post_tpiXSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type post_tpiX_marginal(post_tpiX_marginalSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type post_tmuY(post_tmuYSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type post_tsigma2Y(post_tsigma2YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type post_tpiY(post_tpiYSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type tauY(tauYSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type tauZ(tauZSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type log_pi(log_piSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type log_minus_pi(log_minus_piSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type nuXmat(nuXmatSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type nuYmat(nuYmatSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type post_a0(post_a0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type post_b0(post_b0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type post_a1(post_a1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type post_b1(post_b1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type post_a2x(post_a2xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type post_b2x(post_b2xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type post_a2y(post_a2ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type post_b2y(post_b2ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type meanFactors(meanFactorsSEXP);
    Rcpp::traits::input_parameter< const int >::type seed0(seed0SEXP);
    Rcpp::traits::input_parameter< const double >::type robust_eps(robust_epsSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const double >::type L(LSEXP);
    Rcpp::traits::input_parameter< const double >::type L2(L2SEXP);
    rcpp_result_gen = Rcpp::wrap(spear_(family, Y, X, Yobs, Xobs, Z, nclasses, functional_path, pattern_samples, pattern_features, weights, weights0, weights_case, num_factors, warm_up, max_iter, thres_elbo, thres_count, thres_factor, a0, b0, a1, b1, a2, b2, lower, print_out, interceptsX, interceptsY, post_mu, post_sigma2, post_pi, post_tmuX, post_tsigma2X, post_tpiX, post_tpiX_marginal, post_tmuY, post_tsigma2Y, post_tpiY, tauY, tauZ, log_pi, log_minus_pi, nuXmat, nuYmat, post_a0, post_b0, post_a1, post_b1, post_a2x, post_b2x, post_a2y, post_b2y, meanFactors, seed0, robust_eps, alpha0, L, L2));
    return rcpp_result_gen;
END_RCPP
}
