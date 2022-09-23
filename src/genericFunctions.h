// mixture.h
// =============================================================================
// include guard
#ifndef GENFUN_H
#define GENFUN_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>

using namespace arma ;

// =============================================================================
// a set of functions used in a few different places

//' @title Propose new non-negative value
//' @description Propose new non-negative for sampling.
//' @param x Current value to be proposed
//' @param window The proposal window
//' @return new double
double proposeNewNonNegativeValue(double x, double window, 
  bool use_log_norm = false,
  double tolerance = 1e-8
);

//' @title The Inverse Gamma Distribution
//' @description Random generation from the inverse Gamma distribution.
//' @param shape Shape parameter.
//' @param rate Rate parameter.
//' @return Sample from invGamma(shape, rate).
double rInvGamma(double shape, double rate);

//' @title The Inverse Gamma Distribution
//' @description Random generation from the inverse Gamma distribution.
//' @param N Number of samples to draw.
//' @param shape Shape parameter.
//' @param rate Rate parameter.
//' @return Sample from invGamma(shape, rate).
arma::vec rInvGamma(uword N, double shape, double rate);

//' @title The Gamma Distribution
//' @description Random generation from the Gamma distribution.
//' @param shape Shape parameter.
//' @param rate Rate parameter.
//' @return Sample from Gamma(shape, rate).
double rGamma(double shape, double rate);

//' @title The Gamma Distribution
//' @description Random generation from the Gamma distribution.
//' @param N Number of samples to draw.
//' @param shape Shape parameter.
//' @param rate Rate parameter.
//' @return N samples from Gamma(shape, rate).
arma::vec rGamma(uword N, double shape, double rate);

// //' @title The Inverse Gamma Distribution
// //' @description Random generation from the inverse Gamma distribution.
// //' @param shape Shape parameter.
// //' @param scale Scale parameter.
// //' @return Sample from invGamma(shape, scale).
// double rInvGamma(double shape, double scale);

//' @title The Half-Cauchy Distribution
//' @description Random generation from the Half-Cauchy distribution.
//' See https://en.wikipedia.org/wiki/Cauchy_distribution#Related_distributions
//' @param mu Location parameter.
//' @param scale Scale parameter.
//' @return Sample from HalfCauchy(mu, scale).
double rHalfCauchy(double mu, double scale) ;

//' @title The Beta Distribution
//' @description Random generation from the Beta distribution.
//' See https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions.
//' Samples from a Beta distribution based using two independent gamma
//' distributions.
//' @param a Shape parameter.
//' @param b Shape parameter.
//' @return Sample from Beta(a, b).
double rBeta(double a, double b);

//' @title The Beta Distribution
//' @description Random generation from the Beta distribution.
//' See https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions.
//' Samples from a Beta distribution based using two independent gamma
//' distributions.
//' @param n The number of samples to draw.
//' @param a Shape parameter.
//' @param b Shape parameter.
//' @return Sample from Beta(a, b).
arma::vec rBeta(arma::uword n, double a, double b);

//' @title Metropolis acceptance step
//' @description Given a probaility, randomly accepts by sampling from a uniform 
//' distribution.
//' @param acceptance_prob Double between 0 and 1.
//' @return Boolean indicating acceptance.
bool metropolisAcceptanceStep(double acceptance_prob);

//' @title Squared exponential function
//' @description The squared exponential function as used in a covariance kernel.
//' @param amplitude The amplitude parameter (double)
//' @param length The length parameter (double)
//' @param i Time point (unsigned integer)
//' @param j Time point (unsigned integer)
//' @return Boolean indicating acceptance.
double squaredExponentialFunction(double amplitude, double length, int i, int j);

//' @title Double approximately equal
//' @description Compare two doubles in a way that makes sense.
//' @param x first double considered
//' @param y double compared to x
//' @param precision double of the tolerance of disagreement between x and y.
//' @return bool indicating if the absolute difference between x and y is less 
//' than precision.
bool doubleApproxEqual(double x, double y, double precision = 0.000002);

//' @title Sample mean
//' @description calculate the sample mean of a matrix X.
//' @param X Matrix
//' @return Vector of the column means of X.
arma::vec sampleMean(arma::mat X);

//' @title Calculate sample covariance
//' @description Returns the unnormalised sample covariance. Required as
//' arma::cov() does not work for singletons.
//' @param data Data in matrix format
//' @param sample_mean Sample mean for data
//' @param n The number of samples in data
//' @param n_col The number of columns in data
//' @return One of the parameters required to calculate the posterior of the
//'  Multivariate normal with unknown mean and covariance (the unnormalised
//'  sample covariance).
arma::mat calcSampleCov(
    arma::mat data,
    arma::vec sample_mean,
    arma::uword N,
    arma::uword P
);

//' @title Round matrix
//' @description Round a matrix to n_places decimal places.
//' @param X Matrix
//' @param n_places Integer - number of decimal places to round to
//' @return Matrix X round to n_places decimal places.
arma::mat roundMatrix(arma::mat X, int n_places = 0);

//' @title Choose
//' @description N choose K for binomial coefficient
//' @param n unsigned int (greater than k)
//' @param k unsigned int 
//' @return n choose k
int choose(arma::uword n, arma::uword k);

//' @title Log Choose
//' @description Log transform of N choose K for binomial coefficient
//' @param n unsigned int (greater than k)
//' @param k unsigned int 
//' @return n choose k
double logChoose(double n, double k);

#endif /* GENFUN_H */
