// gaussian.cpp
// =============================================================================
// included dependencies
# include "logLikelihoods.h"
# include "gaussian.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// gaussian class

gaussian::gaussian(arma::uword _K, arma::uvec _labels, arma::mat _X) :
  density(_K, _labels, _X)
{
  mu.zeros(P, K);
  variances.ones(P, K);
  precisions.ones(P, K);
  log_precisions.zeros(P, K);

  // Mean and variance per feature
  n_param = 2 * P;

  kappa = 0.01;
  nu = 3.0;

  // These use only the observed entries
  empiricalBayesHyperparameters();

  identifyMissingValues();
  initializeMissingValues();
};

arma::vec gaussian::empiricalMean() {
  return sampleMeanRobust(X);
};

arma::vec gaussian::empiricalScaleVector() {
  const arma::mat global_cov = computeCovarianceRobust(X);
  const double scale_entry = (arma::accu(global_cov.diag()) / (double) P)
    / std::pow((double) K, 2.0 / (double) P);
  return arma::vec(P, arma::fill::value(scale_entry));
};

void gaussian::empiricalBayesHyperparameters() {
  xi = empiricalMean();
  scale = empiricalScaleVector();
}

Rcpp::List gaussian::hyperparameterList() const {
  return Rcpp::List::create(
    Rcpp::Named("xi") = xi,
    Rcpp::Named("scale") = scale,
    Rcpp::Named("kappa") = kappa,
    Rcpp::Named("nu") = nu
  );
}

void gaussian::setKthVariance(uword k, uword p, double variance) {
  variances(p, k) = variance;
  precisions(p, k) = 1.0 / variance;
  log_precisions(p, k) = -std::log(variance);
}

void gaussian::sampleVariancePrior() {
  for(uword k = 0; k < K; k++){
    for(uword p = 0; p < P; p++) {
      setKthVariance(k, p, 1.0 / rGamma(0.5 * nu, 0.5 * scale(p)));
    }
  }
};

void gaussian::sampleMuPrior() {
  for(uword k = 0; k < K; k++){
    for(uword p = 0; p < P; p++) {
      mu(p, k) = xi(p) + std::sqrt(variances(p, k) / kappa) * randn();
    }
  }
};

void gaussian::sampleFromPriors() {
  sampleVariancePrior();
  sampleMuPrior();
};

// The log likelihood of the observed entries of an item in a component
double gaussian::logLikelihood(arma::uword n, arma::uword k) {
  const arma::uvec& obs_idx = observed_indices(n);
  double ll = 0.0;
  for(uword i = 0; i < obs_idx.n_elem; i++) {
    const uword p = obs_idx(i);
    const double diff = X(n, p) - mu(p, k);
    ll -= 0.5 * (std::log(2.0 * M_PI) - log_precisions(p, k) + precisions(p, k) * diff * diff);
  }
  return ll;
}

arma::vec gaussian::itemLogLikelihood(arma::uword n) {
  arma::vec ll(K);
  for(uword k = 0; k < K; k++) {
    ll(k) = logLikelihood(n, k);
  }
  return ll;
}

void gaussian::sampleKthComponentParameters(
    uword k,
    const umat& members,
    const uvec& non_outliers
) {

  const uvec rel_inds = find((members.col(k) == 1) && (non_outliers == 1));
  const uword n_k = rel_inds.n_elem;

  if(n_k > 0){
    const arma::mat component_data = X.rows( rel_inds ) ;
    const arma::vec sample_mean = mean(component_data, 0).t();
    const double kappa_n = kappa + (double) n_k, nu_n = nu + (double) n_k;
    const arma::vec mu_n = (xi * kappa + (double) n_k * sample_mean) / kappa_n;

    for(uword p = 0; p < P; p++) {
      const double sum_sq = arma::accu(arma::square(component_data.col(p) - sample_mean(p)));
      const double scale_np = scale(p) + sum_sq
        + ((double) n_k * kappa / kappa_n) * std::pow(sample_mean(p) - xi(p), 2.0);

      // sigma^2 | . ~ InvGamma(nu_n / 2, scale_np / 2); mu | sigma^2, . ~ N(mu_n, sigma^2 / kappa_n)
      setKthVariance(k, p, 1.0 / rGamma(0.5 * nu_n, 0.5 * scale_np));
      mu(p, k) = mu_n(p) + std::sqrt(variances(p, k) / kappa_n) * randn();
    }
  } else{
    // Empty components are drawn from the prior
    for(uword p = 0; p < P; p++) {
      setKthVariance(k, p, 1.0 / rGamma(0.5 * nu, 0.5 * scale(p)));
      mu(p, k) = xi(p) + std::sqrt(variances(p, k) / kappa) * randn();
    }
  }
};

void gaussian::sampleMissingForObservation(arma::uword n) {
  const uword k = labels(n);
  const arma::uvec& miss_idx = missing_indices(n);
  for(uword i = 0; i < miss_idx.n_elem; i++) {
    const uword p = miss_idx(i);
    X(n, p) = mu(p, k) + std::sqrt(variances(p, k)) * randn();
  }
}

void gaussian::swapComponents(uword k, uword kprime) {
  mu.swap_cols(k, kprime);
  variances.swap_cols(k, kprime);
  precisions.swap_cols(k, kprime);
  log_precisions.swap_cols(k, kprime);
}

arma::vec gaussian::parameters() const {
  return join_cols(vectorise(mu), vectorise(variances));
}

void gaussian::setParameters(const arma::vec& theta) {
  if(theta.n_elem != 2 * P * K) {
    Rcpp::stop("gaussian: parameter vector has the wrong length.");
  }
  mu = reshape(theta.subvec(0, P * K - 1), P, K);
  variances = reshape(theta.subvec(P * K, 2 * P * K - 1), P, K);
  precisions = 1.0 / variances;
  log_precisions = -log(variances);
}

arma::vec gaussian::simulate(arma::uword k) const {
  return mu.col(k) + sqrt(variances.col(k)) % randn<arma::vec>(P);
}
