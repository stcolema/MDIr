// mvt.cpp
// =============================================================================
# include <RcppArmadillo.h>
# include "mvt.h"

using namespace arma ;

mvt::mvt(
    arma::uvec _fixed, 
    arma::mat _X,
    const arma::field<arma::uvec>* miss_idx,
    const arma::field<arma::uvec>* obs_idx
) : outlierComponent(_fixed, _X, miss_idx, obs_idx) {
  global_mean = sampleMeanRobust(X);
  global_cov = 0.5 * computeCovarianceRobust(X);
  calculateAllLogLikelihoods();
};

// Marginal t density of the observed entries (a sub-vector of a multivariate t
// is multivariate t with the corresponding sub-location and sub-scale)
double mvt::calculateItemLogLikelihood(arma::uword n) {
  const arma::uvec& obs_idx = (*observed_indices_ref)(n);
  if(obs_idx.n_elem == 0) {
    return 0.0;
  }
  arma::vec x_obs = X.row(n).t();
  x_obs = x_obs.elem(obs_idx);
  return mvtLogLikelihood(
    x_obs, 
    global_mean.elem(obs_idx), 
    global_cov.submat(obs_idx, obs_idx), 
    df
  );
}

double mvt::completeLogDensity(const arma::vec& x) const {
  return mvtLogLikelihood(x, global_mean, global_cov, df);
}

// If x ~ t_df(mu, S) then x_m | x_o ~ t_{df + p_o}(mu_m|o, S_m|o (df + d_o) / (df + p_o)) 
// with d_o the squared Mahalanobis distance of x_o.
arma::vec mvt::sampleMissingValues(arma::uword n) const {
  const arma::uvec& miss_idx = (*missing_indices_ref)(n);
  const arma::uvec& obs_idx = (*observed_indices_ref)(n);
  
  arma::vec x_obs;
  if(obs_idx.n_elem > 0) {
    x_obs = X.row(n).t();
    x_obs = x_obs.elem(obs_idx);
  }
  
  arma::vec cond_mean;
  arma::mat cond_cov;
  double mahalanobis = 0.0;
  conditionalMVN(global_mean, global_cov, obs_idx, miss_idx, x_obs, 
    cond_mean, cond_cov, mahalanobis);
  
  const double df_cond = df + (double) obs_idx.n_elem;
  cond_cov *= (df + mahalanobis) / df_cond;
  
  const double mixing = rGamma(0.5 * df_cond, 0.5 * df_cond);
  const arma::vec z = rmvnormChol(arma::zeros<arma::vec>(miss_idx.n_elem), cond_cov);
  return cond_mean + z / std::sqrt(mixing);
}

arma::vec mvt::simulate() const {
  const double mixing = rGamma(0.5 * df, 0.5 * df);
  const arma::vec z = rmvnormChol(arma::zeros<arma::vec>(P), global_cov);
  return global_mean + z / std::sqrt(mixing);
}
