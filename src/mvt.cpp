// mvt.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "mvt.h"

using namespace arma ;

// =============================================================================
// mvt class of outlier component


// Parametrised class
mvt::mvt(arma::uvec _fixed, arma::mat _X) : outlierComponent(_fixed, _X) {
  
  // Doubles for leading coefficient of outlier distribution likelihood
  double lgamma_df_p = 0.0, lgamma_df = 0.0, log_pi_df = 0.0;
  
  // for use in the outlier distribution
  global_cov = findInvertibleGlobalCov();
  global_mean = sampleMean(X);

  // Functions of the covariance relevant to the likelihood
  global_log_det = log_det(global_cov).real();
  global_cov_inv = inv(global_cov);
  
  // Components of the t-distribution likelihood that do not change
  lgamma_df_p = std::lgamma(0.5 * (df + (double) P));
  lgamma_df = std::lgamma(0.5 * df);
  log_pi_df = (0.5 * (double) P) * std::log(df * M_PI);
  
  // Constant in t likelihood
  t_likelihood_const = lgamma_df_p - lgamma_df - log_pi_df - 0.5 * global_log_det;

  // Calculate the log likelihood of each item within the outlier component
  calculateAllLogLikelihoods();
  
};

double mvt::calculateItemLogLikelihood(arma::vec x) {
  
  return mvtLogLikelihood(x, global_mean, global_cov, df);

  // double exponent = 0.0, ll = 0.0;
  // 
  // vec diff_with_mean = x - global_mean;
  // 
  // exponent = as_scalar( diff_with_mean.t() * global_cov_inv * diff_with_mean );
  // 
  // // The T likelihood constant is calculated a member of the TAGM class
  // ll = t_likelihood_const
  //   - 0.5 * (df + (double) P) * std::log(1.0 + (1.0 / df) * exponent);
  // 
  // return ll;
};

arma::mat mvt::findInvertibleGlobalCov(double threshold) {
  
  bool not_invertible = false;
  
  mat small_identity(P, P), global_cov(P, P);
  small_identity.zeros(), global_cov.zeros();
  
  small_identity.eye(P, P);
  small_identity *= 1e-10;
  
  // for use in the outlier distribution
  global_cov = 0.5 * arma::cov(X);
  
  // Do we need to add a very little to the diagonal to ensure we can inverse 
  // the dataset covariance matrix?
  // uword count_here = 0;
  
  vec eigval = eig_sym( global_cov );
  
  not_invertible = min(eigval) < threshold;
  
  // If our covariance matrix is poorly behaved (i.e. non-invertible), add a 
  // small constant to the diagonal entries
  if(not_invertible) {
    global_cov = 0.5 * arma::cov(X) + small_identity;
  }
  
  return global_cov;
};
