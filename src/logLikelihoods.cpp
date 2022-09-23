// logLikelihoods.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "logLikelihoods.h"

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

double wishartLogLikelihood(arma::mat X, arma::mat V, double n, arma::uword P){
  return 0.5*(
    (n - P - 1) * arma::log_det(X).real() 
    - trace(arma::inv_sympd(V) * X) 
    - n * arma::log_det(V).real()
  );
}

double invWishartLogLikelihood(arma::mat X, arma::mat Psi, double nu, arma::uword P) {
  return -0.5 * (
    nu * arma::log_det(Psi).real()
    + (nu + P + 1) * arma::log_det(X).real()
    + arma::trace( Psi * arma::inv_sympd(X) ) 
  );
}

double mvtLogLikelihood(arma::vec x, arma::vec mu, arma::mat Sigma, double nu) {
    
    double P = (double) x.n_rows, exponent = 0.0, ll = 0.0;
    arma::vec mean_diff = x - mu;
    
    exponent = arma::as_scalar(
      mean_diff.t()
      * inv(Sigma)
      * mean_diff
    );

    ll = lgamma(0.5 * (nu + P)) 
      - lgamma(0.5 * nu) 
      - 0.5 * P * log(nu * M_PI)
      - 0.5 * arma::log_det(Sigma).real()
      - ((nu + (double) P) / 2.0) * std::log(1.0 + (1.0 / nu) * exponent);
    
    return ll;
}


double gaussianLogLikelihood(arma::vec x, arma::vec mu, arma::vec std_dev) {
  int P = x.n_rows;
  double ll = 0.0, ll_p = 0.0;
  for(int p = 0; p < P; p++) {
    ll_p = -0.5 * (
      log(2.0 * M_PI) 
      + log(std_dev(p)) 
      + std::pow(x(p) - mu(p), 2.0) / std_dev(p)
    );
    ll += ll_p;
  }
  return ll;
};

double pNorm(arma::vec x, arma::vec mu, arma::mat Sigma, bool is_sympd) {
  // bool cov_is_sympd = Sigma.is_sympd();
  int P = x.n_rows;
  double out = 0.0;
  arma::vec mean_diff = x - mu;
  if(is_sympd) {
    out = -0.5 * (
      (double) P * std::log(2.0 * M_PI) 
      + arma::log_det_sympd(Sigma) 
      + arma::as_scalar(mean_diff.t() * arma::inv_sympd(Sigma) * mean_diff)
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



//' title The Half-Cauchy Distribution
//' description Calculates the pdf of the Half-Cauchy distribution for value x.
//' See https://en.wikipedia.org/wiki/Cauchy_distribution#Related_distributions
//' param x Value to calculate the probability density of.
//' param mu Location parameter.
//' param scale Scale parameter.
//' return Sample from HalfCauchy(mu, scale).
double pHalfCauchy(double x, double mu, double scale, bool logValue) {
  double denom = 0.0;
  
  if(x < mu) {
    Rcpp::Rcerr << "\nIn Half-Cauchy p.d.f, the considered value is less than the threshold.";
    return 0;
  }
  denom = 1.0 + std::pow((x - mu) / scale, 2.0);
  if(logValue) {
    // denom = 2.0 * std::log((x - mu) / scale);
    return log(2.0) - log(M_PI) - log(scale) - log(denom);
  } else {
    return 2.0 / (M_PI * scale * denom);
  }
};
