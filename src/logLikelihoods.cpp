// logLikelihoods.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "logLikelihoods.h"
# include "genericFunctions.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std ;

// =============================================================================

double gammaLogLikelihood(double x, double shape, double rate){
  return shape * log(rate) - lgamma(shape) + (shape - 1) * log(x) - rate * x;
};

double gammaLogLikelihood(arma::vec x, double shape, double rate){
  double out = 0.0;
  for (auto & element : x) {
    out += gammaLogLikelihood(element, shape, rate);
  }
  return out;
};

double invGammaLogLikelihood(double x, double shape, double scale) {
  return shape * log(scale) - lgamma(shape) + (-shape - 1) * log(x) - scale / x;
};

namespace {

// Log of the multivariate gamma function
double logMultivariateGamma(double a, arma::uword P) {
  double out = 0.25 * (double) P * ((double) P - 1.0) * std::log(M_PI);
  for(arma::uword j = 1; j <= P; j++) {
    out += std::lgamma(a + 0.5 * (1.0 - (double) j));
  }
  return out;
}

double logDetSympd(const arma::mat& S) {
  return arma::log_det_sympd(S);
}

}

double wishartLogLikelihood(arma::mat X, arma::mat V, double n, arma::uword P){
  return 0.5 * (
    (n - (double) P - 1.0) * logDetSympd(X)
    - arma::trace(arma::solve(V, X, arma::solve_opts::likely_sympd)) 
    - n * logDetSympd(V)
    - n * (double) P * std::log(2.0)
  ) - logMultivariateGamma(0.5 * n, P);
}

double invWishartLogLikelihood(arma::mat X, arma::mat Psi, double nu, arma::uword P) {
  return 0.5 * (
    nu * logDetSympd(Psi)
    - (nu + (double) P + 1.0) * logDetSympd(X)
    - arma::trace(arma::solve(X, Psi, arma::solve_opts::likely_sympd)) 
    - nu * (double) P * std::log(2.0)
  ) - logMultivariateGamma(0.5 * nu, P);
}

double mvtLogLikelihood(arma::vec x, arma::vec mu, arma::mat Sigma, double nu) {
  const double P = (double) x.n_rows;
  const arma::mat Lower = cholLowerRobust(Sigma);
  const arma::vec z = arma::solve(arma::trimatl(Lower), x - mu);
  const double exponent = arma::dot(z, z);
  
  return std::lgamma(0.5 * (nu + P)) 
    - std::lgamma(0.5 * nu) 
    - 0.5 * P * std::log(nu * M_PI)
    - arma::accu(arma::log(Lower.diag()))
    - 0.5 * (nu + P) * std::log1p(exponent / nu);
}

double pNorm(arma::vec x, arma::vec mu, arma::mat Sigma, bool is_sympd) {
  int P = x.n_rows;
  double out = 0.0;
  arma::vec mean_diff = x - mu;
  if(is_sympd) {
    const arma::mat Lower = cholLowerRobust(Sigma);
    const arma::vec z = arma::solve(arma::trimatl(Lower), mean_diff);
    out = -0.5 * (
      (double) P * std::log(2.0 * M_PI) 
      + 2.0 * arma::accu(arma::log(Lower.diag()))
      + arma::dot(z, z)
    );
  } else {
    out = -0.5 * (
      (double) P * std::log(2.0 * M_PI) 
      + arma::log_det(Sigma).real() 
      + arma::as_scalar(mean_diff.t() * arma::inv(Sigma) * mean_diff)
    );
  }
  return out;
}

double pNorm(double x, double mu, double sigma_2) {
  return (-0.5 * (log(2.0 * M_PI) + log(sigma_2) + pow(x - mu, 2.0) / sigma_2));
}

