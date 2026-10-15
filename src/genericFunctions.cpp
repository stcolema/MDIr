
# include "genericFunctions.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

//' title Propose new non-negative value
//' description Propose new non-negative for sampling.
//' param x Current value to be proposed
//' param window The proposal window
//' return new double
double proposeNewNonNegativeValue(
    double x, 
    double window, 
    bool use_log_norm,
    double tolerance
  ) {
  bool value_below_tolerance = false;
  double proposed_value = 0.0;
  if(use_log_norm) {
    proposed_value = std::exp(std::log(x) + randn() * window);
  } else {
    proposed_value = rGamma(x * window, window);
  }
  
  // If the value is too small (normally close to 0 or negative somehow)
  value_below_tolerance = (proposed_value < tolerance);
  if(value_below_tolerance) {
    proposed_value = proposeNewNonNegativeValue(x, window, use_log_norm, tolerance);
  }
  
  return proposed_value;
};

//' title The Inverse Gamma Distribution
//' description Random generation from the inverse Gamma distribution.
//' param shape Shape parameter.
//' param rate Rate parameter.
//' return Sample from invGamma(shape, rate).
double rInvGamma(double shape, double rate) {
  double x = arma::randg( distr_param(shape, 1.0 / rate) );
  return (1 / x);
};

//' title The Inverse Gamma Distribution
//' description Random generation from the inverse Gamma distribution.
//' param N Number of samples to draw.
//' param shape Shape parameter.
//' param rate Rate parameter.
//' return Sample from invGamma(shape, rate).
arma::vec rInvGamma(uword N, double shape, double rate) {
  vec x = arma::randg(N, distr_param(shape, 1.0 / rate) );
  return (1 / x);
};

//' title The Gamma Distribution
//' description Random generation from the Gamma distribution.
//' param shape Shape parameter.
//' param rate Rate parameter.
//' return Sample from Gamma(shape, rate).
double rGamma(double shape, double rate) {
  // Gamma draws with small shape can underflow to exactly zero, which breaks
  // any subsequent log(). Floor at the smallest positive normal double.
  double x = arma::randg( distr_param(shape, 1.0 / rate) );
  return std::max(x, std::numeric_limits<double>::min());
};

//' title The Gamma Distribution
//' description Random generation from the Gamma distribution.
//' param N Number of samples to draw.
//' param shape Shape parameter.
//' param rate Rate parameter.
//' return N samples from Gamma(shape, rate).
arma::vec rGamma(uword N, double shape, double rate) {
  return arma::randg(N, distr_param(shape, 1.0 / rate) );
};


//' title The Half-Cauchy Distribution
//' description Random generation from the Half-Cauchy distribution.
//' See https://en.wikipedia.org/wiki/Cauchy_distribution#Related_distributions
//' param mu Location parameter.
//' param scale Scale parameter.
//' return Sample from HalfCauchy(mu, scale).
double rHalfCauchy(double mu, double scale) {
  double x = 0.0, y = 0.0;
  x = randn();
  while(x <= 0.0) {
    x = randn();
  }
  y = rInvGamma(0.5, 0.5 * std::pow(scale, 2.0));
  return mu + x * std::sqrt(y);
};

//' title The Half-Cauchy Distribution
//' description Random generation from the Half-Cauchy distribution.
//' See https://en.wikipedia.org/wiki/Cauchy_distribution#Related_distributions
//' param N The number of samples to draw
//' param mu Location parameter.
//' param scale Scale parameter.
//' return Sample from HalfCauchy(mu, scale).
arma::vec rHalfCauchy(uword N, arma::vec mu, double scale) {
  vec x(N), y(N);
  x = arma::abs(arma::randn(N));
  y = rInvGamma(N, 0.5, 0.5 * std::pow(scale, 2.0));
  return mu + x * arma::sqrt(y);
};

//' title The Beta Distribution
//' description Random generation from the Beta distribution.
//' See https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions.
//' Samples from a Beta distribution based using two independent gamma
//' distributions.
//' param a Shape parameter.
//' param b Shape parameter.
//' return Sample from Beta(a, b).
double rBeta(double a, double b) { // double theta = 1.0) {
  double X = arma::randg( arma::distr_param(a, 1.0) );
  double Y = arma::randg( arma::distr_param(b, 1.0) );
  double beta = X / (double)(X + Y);
  return(beta);
};

//' title The Beta Distribution
//' description Random generation from the Beta distribution.
//' See https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions.
//' Samples from a Beta distribution based using two independent gamma
//' distributions.
//' param n The number of samples to draw.
//' param a Shape parameter.
//' param b Shape parameter.
//' return Sample from Beta(a, b).
arma::vec rBeta(arma::uword n, double a, double b) {
  arma::vec X = arma::randg(n, arma::distr_param(a, 1.0) );
  arma::vec Y = arma::randg(n, arma::distr_param(b, 1.0) );
  arma::vec beta = X / (X + Y);
  return(beta);
};


double logSumExp(const arma::vec& x) {
  const double m = x.max();
  if(!std::isfinite(m)) {
    return m;
  }
  return m + std::log(arma::accu(arma::exp(x - m)));
}

arma::uword sampleCategorical(const arma::vec& probs) {
  const double u = arma::randu();
  double cumulative = 0.0;
  for(arma::uword i = 0; i + 1 < probs.n_elem; i++) {
    cumulative += probs(i);
    if(u < cumulative) {
      return i;
    }
  }
  return probs.n_elem - 1;
}

//' title Metropolis acceptance step
//' description Given a probaility, randomly accepts by sampling from a uniform 
//' distribution.
//' param acceptance_prob Double between 0 and 1.
//' return Boolean indicating acceptance.
bool metropolisAcceptanceStep(double acceptance_prob) {
  double u = arma::randu();
  return (u < acceptance_prob);
};

//' title Squared exponential function
//' description The squared exponential function as used in a covariance kernel.
//' param amplitude The amplitude parameter (double)
//' param length The length parameter (double)
//' param i Time point (unsigned integer)
//' param j Time point (unsigned integer)
//' return Squared exponential metric of (i, j)
double squaredExponentialFunction(double amplitude, double length, int i, int j) {
  // if(i > j) {
  //   return amplitude * std::exp(- std::pow(i - j, 2.0) / length);
  // } 
  return amplitude * std::exp(- std::pow((double) (j - i), 2.0) / (2.0 * length));
};

bool doubleApproxEqual(double x, double y, double precision) {
  return std::abs(x - y) < precision;
};


// title Sample mean
// description calculate the sample mean of a matrix X.
// param X Matrix
// return Vector of the column means of X.
arma::vec sampleMean(arma::mat X) {
  arma::mat mu_t = arma::mean(X);
  return mu_t.row(0).t();
};

// Compute mean robustly with missing values
arma::vec sampleMeanRobust(const arma::mat& X) {
  arma::uword P = X.n_cols;
  arma::vec means(P);
  
  for(arma::uword p = 0; p < P; p++) {
    arma::vec col_data = X.col(p);
    arma::uvec finite_indices = arma::find_finite(col_data);
    
    if(finite_indices.n_elem > 0) {
      means(p) = arma::mean(col_data.elem(finite_indices));
    } else {
      means(p) = 0.0; // Default if all missing
    }
  }
  
  return means;
}

// Covariance from pairwise-complete observations. This is used only to set
// data-driven hyperparameters, so an approximation that tolerates any pattern
// of missingness is what is wanted. The pairwise estimate need not be positive
// semi-definite, so negative eigenvalues are lifted to a small positive floor.
// Columns with fewer than two finite values receive unit variance.
arma::mat computeCovarianceRobust(const arma::mat& X) {
  const arma::uword P = X.n_cols;
  arma::mat S(P, P, arma::fill::zeros);
  arma::vec means = sampleMeanRobust(X);
  
  for(arma::uword p = 0; p < P; p++) {
    for(arma::uword q = p; q < P; q++) {
      double total = 0.0;
      arma::uword n_pair = 0;
      for(arma::uword n = 0; n < X.n_rows; n++) {
        if(std::isfinite(X(n, p)) && std::isfinite(X(n, q))) {
          total += (X(n, p) - means(p)) * (X(n, q) - means(q));
          n_pair++;
        }
      }
      if(n_pair > 1) {
        S(p, q) = total / (double) (n_pair - 1);
      } else if(p == q) {
        S(p, q) = 1.0;
      }
      S(q, p) = S(p, q);
    }
  }
  
  arma::vec eigval;
  arma::mat eigvec;
  if(arma::eig_sym(eigval, eigvec, S)) {
    const double floor_value = std::max(1e-8, 1e-6 * eigval.max());
    if(eigval.min() < floor_value) {
      eigval = arma::clamp(eigval, floor_value, arma::datum::inf);
      S = eigvec * arma::diagmat(eigval) * eigvec.t();
      S = 0.5 * (S + S.t());
    }
  }
  return S;
}

// Lower Cholesky factor of a covariance matrix. If the matrix is numerically
// indefinite a multiple of the identity is added, growing geometrically, until
// the factorisation succeeds.
arma::mat cholLowerRobust(const arma::mat& S) {
  arma::mat A = 0.5 * (S + S.t()), Lower;
  if(arma::chol(Lower, A, "lower")) {
    return Lower;
  }
  double jitter = 1e-10 * std::max(arma::trace(A) / (double) A.n_rows, 1e-12);
  for(int attempt = 0; attempt < 12; attempt++) {
    arma::mat B = A;
    B.diag() += jitter;
    if(arma::chol(Lower, B, "lower")) {
      return Lower;
    }
    jitter *= 10.0;
  }
  Rcpp::stop("Covariance matrix is not positive definite and could not be repaired.");
  return Lower;
}

// Draw from N(mean, cov) via Cholesky.
arma::vec rmvnormChol(const arma::vec& mean, const arma::mat& cov) {
  arma::mat Lower = cholLowerRobust(cov);
  return mean + Lower * arma::randn<arma::vec>(mean.n_elem);
}

// Conditional distribution of the entries `miss` of a N(mu, Sigma) vector given 
// the entries `obs`. Also returns the squared Mahalanobis distance of x_obs 
// from mu_obs (needed for the multivariate t conditional).
void conditionalMVN(
    const arma::vec& mu,
    const arma::mat& Sigma,
    const arma::uvec& obs,
    const arma::uvec& miss,
    const arma::vec& x_obs,
    arma::vec& cond_mean,
    arma::mat& cond_cov,
    double& mahalanobis_obs
) {
  if(obs.n_elem == 0) {
    cond_mean = mu.elem(miss);
    cond_cov = Sigma.submat(miss, miss);
    mahalanobis_obs = 0.0;
    return;
  }
  arma::mat S_oo = Sigma.submat(obs, obs), S_mo = Sigma.submat(miss, obs);
  arma::mat Lower = cholLowerRobust(S_oo);
  arma::vec diff = x_obs - mu.elem(obs);
  arma::vec z = arma::solve(arma::trimatl(Lower), diff);
  mahalanobis_obs = arma::dot(z, z);
  arma::mat B = arma::solve(arma::trimatl(Lower), S_mo.t());     // L^{-1} S_om
  cond_mean = mu.elem(miss) + B.t() * z;
  cond_cov = Sigma.submat(miss, miss) - B.t() * B;
  cond_cov = 0.5 * (cond_cov + cond_cov.t());
}

// title Calculate sample covariance
// description Returns the unnormalised sample covariance. Required as
// arma::cov() does not work for singletons.
// param data Data in matrix format
// param sample_mean Sample mean for data
// param n The number of samples in data
// param n_col The number of columns in data
// return One of the parameters required to calculate the posterior of the
//  Multivariate normal with uknown mean and covariance (the unnormalised
//  sample covariance).
arma::mat calcSampleCov(const arma::mat& data,
                        const arma::vec& sample_mean,
                        arma::uword N,
                        arma::uword P
) {
  
  mat sample_covariance = zeros<mat>(P, P);
  
  // If n > 0 (as this would crash for empty clusters), and for n = 1 the
  // sample covariance is 0
  if(N > 1){
    // Make explicit copy before centering
    arma::mat centered_data = data;
    centered_data.each_row() -= sample_mean.t();
    sample_covariance = centered_data.t() * centered_data;
    // sample_covariance = ((double) N - 1.0) * cov(data);
  }
  return sample_covariance;
};

arma::mat roundMatrix(arma::mat X, int n_places) {
  double multiplier = std::pow(10, n_places);
  return round(X * multiplier) / multiplier;
}

double logChoose(double n, double k) {
  if(k < 0.0 || k > n) {
    return -arma::datum::inf;
  }
  return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
}
