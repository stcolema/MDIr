// mvn.cpp
// =============================================================================
// included dependencies
# include "logLikelihoods.h"
# include "mvn.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// mvn class

mvn::mvn(arma::uword _K, arma::uvec _labels, arma::mat _X, arma::vec _density_prior) :
  density(_K, _labels, _X, _density_prior)
{
  mu.set_size(P, K);
  mu.zeros();

  cov.set_size(P, P, K);
  cov.zeros();

  cov_log_det = arma::zeros<arma::vec>(K);
  cov_inv.set_size(P, P, K);
  cov_inv.zeros();

  // Mean vector and covariance matrix
  n_param = P * (1 + (P + 1) * 0.5);

  kappa = 0.01;
  nu = P + 2;

  // These use only the observed entries
  empiricalBayesHyperparameters();
  
  // The empirical scale is the prior mean of the pooled scale
  scale_shape = density_prior(0);
  scale_prior_mean = scale.diag();

  identifyMissingValues();
  initializeMissingValues();
};

arma::vec mvn::empiricalMean() {
  return sampleMeanRobust(X);
};

arma::mat mvn::empiricalScaleMatrix() {
  const arma::mat global_cov = computeCovarianceRobust(X);
  const double scale_entry = (arma::accu(global_cov.diag()) / (double) P)
    / std::pow((double) K, 2.0 / (double) P);
  return scale_entry * arma::eye(P, P);
};

void mvn::empiricalBayesHyperparameters() {
  xi = empiricalMean();
  scale = empiricalScaleMatrix();
}

Rcpp::List mvn::hyperparameterList() const {
  return Rcpp::List::create(
    Rcpp::Named("xi") = xi,
    Rcpp::Named("scale") = scale,
    Rcpp::Named("scale_prior_mean") = scale_prior_mean,
    Rcpp::Named("scale_pooling_shape") = scale_shape,
    Rcpp::Named("kappa") = kappa,
    Rcpp::Named("nu") = nu
  );
}

void mvn::sampleCovPrior() {
  for(arma::uword k = 0; k < K; k++){
    cov.slice(k) = arma::iwishrnd(scale, nu);
  }
};

void mvn::sampleMuPrior() {
  for(arma::uword k = 0; k < K; k++){
    mu.col(k) = rmvnormChol(xi, (1.0 / kappa) * cov.slice(k));
  }
};

void mvn::updatePooledHyperparameters(const arma::uvec& occupied) {
  if(scale_shape <= 0.0) {
    return;
  }
  const double shape = scale_shape + 0.5 * nu * (double) occupied.n_elem;
  arma::vec s(P);
  for(uword p = 0; p < P; p++) {
    double rate = scale_shape / scale_prior_mean(p);
    for(uword i = 0; i < occupied.n_elem; i++) {
      rate += 0.5 * cov_inv(p, p, occupied(i));
    }
    s(p) = rGamma(shape, rate);
  }
  scale = arma::diagmat(s);
}

arma::vec mvn::pooledHyperparameters() const {
  return scale.diag();
}

void mvn::setPooledHyperparameters(const arma::vec& pooled) {
  if(pooled.n_elem == P) {
    scale = arma::diagmat(pooled);
  }
}

void mvn::sampleFromPriors() {
  if(scale_shape > 0.0) {
    arma::vec s(P);
    for(uword p = 0; p < P; p++) {
      s(p) = rGamma(scale_shape, scale_shape / scale_prior_mean(p));
    }
    scale = arma::diagmat(s);
  }
  sampleCovPrior();
  sampleMuPrior();
  matrixCombinations();
};

void mvn::matrixCombinations() {
  for(arma::uword k = 0; k < K; k++) {
    cov_inv.slice(k) = arma::inv_sympd(cov.slice(k));
    cov_log_det(k) = arma::log_det_sympd(cov.slice(k));
  }
};

void mvn::sampleKthComponentParameters(
    uword k,
    const umat& members,
    const uvec& non_outliers
  ) {

  // Outliers do not contribute to the component parameters
  const uvec rel_inds = find((members.col(k) == 1) && (non_outliers == 1));
  const arma::uword n_k = rel_inds.n_elem;

  if(n_k > 0){
    const arma::mat component_data = X.rows( rel_inds ) ;
    const arma::vec sample_mean = mean(component_data, 0).t();
    const arma::mat sample_cov = calcSampleCov(component_data, sample_mean, n_k, P);
    const arma::vec dist = sample_mean - xi;

    // Tempered conjugate update: the likelihood enters as L^beta, which is the
    // untempered update with n_k replaced by beta n_k in the counts and
    // beta S for the scatter (the sample mean is unchanged).
    const double n_eff = beta * (double) n_k;
    const double kappa_n = kappa + n_eff;
    const arma::mat scale_n = scale + beta * sample_cov
      + ((kappa * n_eff) / kappa_n) * (dist * dist.t());

    cov.slice(k) = iwishrnd(scale_n, nu + n_eff);

    const arma::vec mu_n = (kappa * xi + n_eff * sample_mean) / kappa_n;
    mu.col(k) = rmvnormChol(mu_n, cov.slice(k) / kappa_n);
  } else{
    // Empty components are drawn from the prior
    cov.slice(k) = iwishrnd(scale, nu);
    mu.col(k) = rmvnormChol(xi, cov.slice(k) / kappa);
  }

  cov_inv.slice(k) = inv_sympd(cov.slice(k));
  cov_log_det(k) = log_det_sympd(cov.slice(k));
};

void mvn::sampleMissingForObservation(arma::uword n) {
  const arma::uvec& miss_idx = missing_indices(n);
  if(miss_idx.n_elem == 0) {
    return;
  }
  const arma::uvec& obs_idx = observed_indices(n);
  const uword k = labels(n);

  arma::vec x_obs, cond_mean;
  arma::mat cond_cov;
  double mahalanobis = 0.0;
  if(obs_idx.n_elem > 0) {
    x_obs = X.row(n).t();
    x_obs = x_obs.elem(obs_idx);
  }
  conditionalMVN(mu.col(k), cov.slice(k), obs_idx, miss_idx, x_obs,
    cond_mean, cond_cov, mahalanobis);

  const arma::vec sampled = rmvnormChol(cond_mean, cond_cov);
  for(uword i = 0; i < miss_idx.n_elem; i++) {
    X(n, miss_idx(i)) = sampled(i);
  }
}

namespace {

// (x - mu)' S (x - mu) for a P x P symmetric S held column-major. The
// difference is written into the scratch buffer d (length P). Plain loops avoid
// the temporaries of the matrix expression, which dominate for the small P
// typical of a component.
inline double quadraticForm(
    const double* x,
    const double* mu,
    const double* S,
    arma::uword P,
    double* d
) {
  for(arma::uword i = 0; i < P; i++) {
    d[i] = x[i] - mu[i];
  }
  double q = 0.0;
  for(arma::uword j = 0; j < P; j++) {
    const double* S_j = S + j * P;
    double col = 0.0;
    for(arma::uword i = 0; i < P; i++) {
      col += S_j[i] * d[i];
    }
    q += d[j] * col;
  }
  return q;
}

}

double mvn::logLikelihood(arma::uword n, arma::uword k) {
  const arma::uvec& obs_idx = observed_indices(n);
  const uword P_obs = obs_idx.n_elem;

  if(P_obs == P) {
    // Complete data: use the cached inverse and log determinant
    const arma::vec x = X.row(n).t();
    arma::vec d(P);
    return -0.5 * (
      (double) P * std::log(2.0 * M_PI)
      + cov_log_det(k)
      + quadraticForm(x.memptr(), mu.colptr(k), cov_inv.slice(k).memptr(), P, d.memptr())
    );
  } else if(P_obs > 0) {
    // Marginal density of the observed entries
    arma::vec x_obs = X.row(n).t();
    x_obs = x_obs.elem(obs_idx);
    const arma::vec mu_k = mu.col(k);
    const arma::vec mu_obs = mu_k.elem(obs_idx);
    const arma::mat cov_obs = cov.slice(k).submat(obs_idx, obs_idx);
    return pNorm(x_obs, mu_obs, cov_obs, true);
  }
  // Nothing observed: the likelihood does not depend on the component
  return 0.0;
}

arma::vec mvn::itemLogLikelihood(arma::uword n) {
  arma::vec ll(K);
  if(observed_indices(n).n_elem == P) {
    // Complete data: read the item once rather than once per component
    const arma::vec x = X.row(n).t();
    arma::vec d(P);
    for(uword k = 0; k < K; k++) {
      ll(k) = -0.5 * (
        (double) P * std::log(2.0 * M_PI)
        + cov_log_det(k)
        + quadraticForm(x.memptr(), mu.colptr(k), cov_inv.slice(k).memptr(), P, d.memptr())
      );
    }
    return ll;
  }
  for(uword k = 0; k < K; k++) {
    ll(k) = logLikelihood(n, k);
  }
  return ll;
}

void mvn::swapComponents(uword k, uword kprime) {
  mu.swap_cols(k, kprime);
  cov.slice(k).swap(cov.slice(kprime));
  cov_inv.slice(k).swap(cov_inv.slice(kprime));
  std::swap(cov_log_det(k), cov_log_det(kprime));
}

arma::vec mvn::parameters() const {
  return join_cols(vectorise(mu), vectorise(cov));
}

void mvn::setParameters(const arma::vec& theta) {
  const uword n_mu = P * K;
  if(theta.n_elem != n_mu + P * P * K) {
    Rcpp::stop("mvn: parameter vector has the wrong length.");
  }
  mu = reshape(theta.subvec(0, n_mu - 1), P, K);
  for(uword k = 0; k < K; k++) {
    cov.slice(k) = reshape(theta.subvec(n_mu + k * P * P, n_mu + (k + 1) * P * P - 1), P, P);
  }
  matrixCombinations();
}

arma::vec mvn::simulate(arma::uword k) const {
  return rmvnormChol(mu.col(k), cov.slice(k));
}


// === Collapsed marginal likelihood ============================================
// With n items, beta-tempered likelihood and the NIW(xi, kappa, scale, nu) prior,
//   log m = -(beta n P / 2) log(pi) + (P / 2) log(kappa / (kappa + beta n))
//           + log Gamma_P((nu + beta n) / 2) - log Gamma_P(nu / 2)
//           + (nu / 2) log|scale| - ((nu + beta n) / 2) log|scale_n|,
//   scale_n = scale + beta S + kappa beta n / (kappa + beta n) (xbar - xi)(xbar - xi)',
// S the scatter about the sample mean. Checked against the identity
// prior x L^beta = posterior x m in verification/tempering/sympy_checks.py.
collapsedStats mvn::emptyStats() const {
  collapsedStats st;
  st.n = 0.0;
  st.s.zeros(P);
  st.S2.zeros(P, P);
  return st;
}

void mvn::addItemToStats(collapsedStats& st, arma::uword n) const {
  const arma::vec d = X.row(n).t() - xi;
  st.n += 1.0;
  st.s += d;
  st.S2 += d * d.t();
}

double mvn::logMarginalLikelihood(const collapsedStats& st, double beta) const {
  if(st.n <= 0.0) {
    return 0.0;
  }
  const double n_eff = beta * st.n, kappa_n = kappa + n_eff, Pd = (double) P;
  const arma::vec dbar = st.s / st.n;
  const arma::mat S = st.S2 - st.s * st.s.t() / st.n;
  const arma::mat scale_n = scale + beta * S + (kappa * n_eff / kappa_n) * (dbar * dbar.t());
  auto logMultivariateGamma = [&](double a) {
    double out = Pd * (Pd - 1.0) / 4.0 * std::log(M_PI);
    for(uword j = 0; j < P; j++) {
      out += std::lgamma(a - 0.5 * (double) j);
    }
    return out;
  };
  double logdet_scale = 0.0, logdet_scale_n = 0.0;
  arma::log_det_sympd(logdet_scale, scale);
  arma::log_det_sympd(logdet_scale_n, 0.5 * (scale_n + scale_n.t()));
  return -0.5 * n_eff * Pd * std::log(M_PI)
    + 0.5 * Pd * (std::log(kappa) - std::log(kappa_n))
    + logMultivariateGamma(0.5 * (nu + n_eff)) - logMultivariateGamma(0.5 * nu)
    + 0.5 * nu * logdet_scale - 0.5 * (nu + n_eff) * logdet_scale_n;
}
