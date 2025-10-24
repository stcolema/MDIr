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

// double mvt::calculateItemLogLikelihood(arma::vec x) {
//   
//   return mvtLogLikelihood(x, global_mean, global_cov, df);
// 
//   // double exponent = 0.0, ll = 0.0;
//   // 
//   // vec diff_with_mean = x - global_mean;
//   // 
//   // exponent = as_scalar( diff_with_mean.t() * global_cov_inv * diff_with_mean );
//   // 
//   // // The T likelihood constant is calculated a member of the TAGM class
//   // ll = t_likelihood_const
//   //   - 0.5 * (df + (double) P) * std::log(1.0 + (1.0 / df) * exponent);
//   // 
//   // return ll;
// };

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

void mvt::initializeMissingValues() {
  if(missing_indices_ref == nullptr) return;
  
  for(uword n = 0; n < N; n++) {
    if((*missing_indices_ref)(n).n_elem > 0) {
      arma::uvec miss_idx = (*missing_indices_ref)(n);
      for(uword idx : miss_idx) {
        X(n, idx) = global_mean(idx);  // Initialize with global mean
      }
    }
  }
  X_t = X.t();
}


void mvt::sampleMissingForObservation(arma::uword n) {
  if(missing_indices_ref == nullptr || (*missing_indices_ref)(n).n_elem == 0) return;
  
  arma::uvec miss_idx = (*missing_indices_ref)(n);
  arma::uvec obs_idx = (*observed_indices_ref)(n);
  
  if(obs_idx.n_elem > 0) {
    // CORRECT: Use .elem() on full vectors
    arma::vec x_n_full = X.row(n).t();
    
    arma::vec x_obs = x_n_full.elem(obs_idx);
    arma::vec mu_obs = global_mean.elem(obs_idx);
    arma::vec mu_miss = global_mean.elem(miss_idx);
    
    arma::mat cov_miss = global_cov.submat(miss_idx, miss_idx);
    arma::mat cov_obs = global_cov.submat(obs_idx, obs_idx);
    arma::mat cov_cross = global_cov.submat(miss_idx, obs_idx);
    
    arma::vec conditional_mean = mu_miss + cov_cross * arma::solve(cov_obs, x_obs - mu_obs);
    arma::mat conditional_cov = cov_miss - cov_cross * arma::solve(cov_obs, cov_cross.t());
    
    arma::vec sampled = arma::mvnrnd(conditional_mean, conditional_cov);
    for(uword i = 0; i < miss_idx.n_elem; i++) {
      X(n, miss_idx(i)) = sampled(i);
    }
  } else {
    // All missing
    arma::vec mu_miss = global_mean.elem(miss_idx);
    arma::mat cov_miss = global_cov.submat(miss_idx, miss_idx);
    arma::vec sampled = arma::mvnrnd(mu_miss, cov_miss);
    for(uword i = 0; i < miss_idx.n_elem; i++) {
      X(n, miss_idx(i)) = sampled(i);
    }
  }
}
  
double mvt::calculateItemLogLikelihood(arma::uword n) {
  if(observed_indices_ref == nullptr) {
    return mvtLogLikelihood(X.row(n).t(), global_mean, global_cov, df);
  }
  
  arma::uvec obs_idx = (*observed_indices_ref)(n);
  
  if(obs_idx.n_elem == P) {
    return mvtLogLikelihood(X.row(n).t(), global_mean, global_cov, df);
  } else if(obs_idx.n_elem > 0) {
    // CORRECT: Use .elem() on full vectors
    arma::vec x_n_full = X.row(n).t();
    arma::vec x_obs = x_n_full.elem(obs_idx);
    arma::vec mu_obs = global_mean.elem(obs_idx);
    arma::mat cov_obs = global_cov.submat(obs_idx, obs_idx);
    
    return mvtLogLikelihood(x_obs, mu_obs, cov_obs, df);
  } else {
    return 0.0;
  }
}