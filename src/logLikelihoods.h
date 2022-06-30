// logLikelihoods.h
// =============================================================================
  // include guard
#ifndef LOGLIKELIHOOD_H
#define LOGLIKELIHOOD_H

// =============================================================================
  // included dependencies
# include <RcppArmadillo.h>

// =============================================================================
  
//' @title Gamma log-likelihood
//' @description Used in calculating model probability in Metropolis-Hastings 
//' algorithm when proposals are from the Gamma distribution.
//' @param x - double; the value to calculate the likelihood of.
//' @param shape - double; the shape of the Gamma distribution.
//' @param rate - double; the rate of the Gamma distribution
//' @return the log-likelihood of x in a Gamma with parameters shape 
//' and rate.
// [[Rcpp::export]]
double gammaLogLikelihood(double x, double shape, double rate);

//' @title Gamma log-likelihood
//' @description The log-likelihood of each element of a vector in a Gamma 
//' distribution parametrised with a shape and rate.
//' @param x - vector; the values to calculate the likelihood of.
//' @param shape - double; the shape of the Gamma distribution.
//' @param rate - double; the rate of the Gamma distribution
//' @return the log-likelihood of x in a Gamma with parameters shape and rate.
double gammaLogLikelihood(arma::vec x, double shape, double rate);

//' @title Inverse gamma log-likelihood
//' @description Used in calculating model probability in Metropolis-Hastings 
//' algorithm when proposals are from the inverse-Gamma distribution.
//' @param x - double; the value to calculate the likelihood of.
//' @param shape - double; the shape of the inverse-Gamma distribution.
//' @param scale - double; the scale of the inverse-Gamma distribution
//' @return the unnormalised log-likelihood of x in a inverse-Gamma with parameters 
//' shape and scale.
// [[Rcpp::export]]
double invGammaLogLikelihood(double x, double shape, double scale);

//' @title Wishart log-likelihood
//' @description Used in calculating model probability in Metropolis-Hastings 
//' algorithm when proposals are from the Wishart distribution.
//' @param X - matrix; the matrix to calculate the likelihood of.
//' @param V - matrix; the scale of the Wishart distribution.
//' @param n - double; the degrees of freedom for the Wishart distribution.
//' @param P - unsigned integer; the dimension of X.
//' @return the unnormalised log-likelihood of X in a Wishart with parameters V 
//' and n.
// [[Rcpp::export]]
double wishartLogLikelihood(arma::mat X, arma::mat V, double n, arma::uword P);

//' @title Inverse-Wishart log-likelihood
//' @description Used in calculating model probability in Metropolis-Hastings 
//' algorithm when proposals are from the Wishart distribution.
//' @param X - matrix; the matrix to calculate the likelihood of.
//' @param Psi - matrix; the scale of the inverse-Wishart distribution.
//' @param nu - double; the degrees of freedom for the inverse-Wishart distribution.
//' @param P - unsigned integer; the dimension of X.
//' @return the unnormalised log-likelihood of X in a inverse-Wishart with parameters Psi 
//' and nu.
// [[Rcpp::export]]
double invWishartLogLikelihood(arma::mat X, arma::mat Psi, double nu, arma::uword P);

//' @title Multivariate t log-likelihood
//' @description The log-likelihood function for a point in the multivariate t 
//' (MVT) distribution.
//' @param x - vector; the sample to calculate the log likelihood of.
//' @param mu - vector; the mean parameter of the MVT distribution.
//' @param Sigma - matrix; the scale matrix of the MVT distribution.
//' @param nu - double; the degrees of freedom for the MVT distribution.
//' @return the normalised log-likelihood of x in a MVT distribution with 
//' parameters mu, Sigma and nu.
// [[Rcpp::export]]
double mvtLogLikelihood(arma::vec x, arma::vec mu, arma::mat Sigma, double nu);

//' @title Multivariate Normal log-likelihood
//' @description The log-likelihood function for a point in the multivariate 
//' Normal (MVN) distribution.
//' @param x - vector; the sample to calculate the log likelihood of.
//' @param mu - vector; the mean parameter of the MVN distribution.
//' @param Sigma - matrix; the covariance matrix of the MVN distribution.
//' @param is_sympd - boolean; is the covariance matrix positive definite (
//' calculations are faster if this is the case).
//' @return the normalised log-likelihood of x in a MVN distribution with 
//' parameters mu, Sigma.
// [[Rcpp::export]]
double pNorm(arma::vec x, arma::vec mu, arma::mat Sigma, bool is_sympd = true);

//' @title Gaussian log-likelihood
//' @description The log-likelihood function for a point in the univariate 
//' Gaussian distribution.
//' @param x - double; the sample to calculate the log likelihood of.
//' @param mu - double; the mean parameter of the Gaussian distribution.
//' @param sigma_2 - double; the standard deviation of the Gaussian distribution.
//' @return the normalised log-likelihood of x in a Gaussian distribution with 
//' parameters mu, sigma_2.
double pNorm(double x, double mu, double sigma_2);

#endif /* LOGLIKELIHOOD_H */