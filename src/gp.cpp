// gp.cpp
// =============================================================================
// included dependencies
# include "logLikelihoods.h"
# include "gp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// gp class

gp::gp(arma::uword _K, arma::uvec _labels, arma::mat _X) : 
  density(_K, _labels, _X) 
{
  amplitude.ones(K);
  length.ones(K);
  noise.ones(K);
  mu.zeros(P, K);
  kernel_sub_block.zeros(P, P, K);
  I_p = eye(P, P);
  
  // -(i - j)^2 / 2, so that K = amplitude * exp(time_diff_mat / length)
  time_diff_mat.zeros(P, P);
  for(uword ii = 0; ii < P; ii++) {
    for(uword jj = ii + 1; jj < P; jj++) {
      time_diff_mat(ii, jj) = - 0.5 * std::pow((double) (jj - ii), 2.0);
      time_diff_mat(jj, ii) = time_diff_mat(ii, jj);
    }
  }
  
  noise_acceptance_count.zeros(K);
  length_acceptance_count.zeros(K);
  amplitude_acceptance_count.zeros(K);
  
  // Mean function plus a noise per component
  n_param = P + 3;
  
  hypers.zeros(3 * K);
  acceptance_count.zeros(3 * K);
  
  identifyMissingValues();
  initializeMissingValues();
};

Rcpp::List gp::hyperparameterList() const {
  return Rcpp::List::create(
    Rcpp::Named("hyper_prior_sd") = hyper_prior_sd,
    Rcpp::Named("noise_prior_sd") = noise_prior_sd,
    Rcpp::Named("kernel_jitter") = kernel_jitter
  );
}

void gp::recordHypers() {
  hypers.subvec(0, K - 1) = amplitude;
  hypers.subvec(K, 2 * K - 1) = length;
  hypers.subvec(2 * K, 3 * K - 1) = noise;
  acceptance_count.subvec(0, K - 1) = amplitude_acceptance_count;
  acceptance_count.subvec(K, 2 * K - 1) = length_acceptance_count;
  acceptance_count.subvec(2 * K, 3 * K - 1) = noise_acceptance_count;
}

// === Priors ==================================================================

// A draw from the log-normal prior restricted to the permitted range
double gp::sampleHyperPrior(double sd) const {
  double x = 0.0;
  do {
    x = std::exp(sd * randn());
  } while(x < hyper_lower || x > hyper_upper);
  return x;
}

void gp::sampleKthComponentHyperParameterPrior(uword k) {
  amplitude(k) = sampleHyperPrior(hyper_prior_sd);
  length(k) = sampleHyperPrior(hyper_prior_sd);
  noise(k) = sampleHyperPrior(noise_prior_sd);
  kernel_sub_block.slice(k) = calculateKthComponentKernelSubBlock(amplitude(k), length(k));
};

void gp::sampleFromPriors() {
  for(uword k = 0; k < K; k++) {
    sampleKthComponentHyperParameterPrior(k);
    mu.col(k) = rmvnormChol(zeros<vec>(P), kernel_sub_block.slice(k));
  }
  recordHypers();
};

// === Kernel ==================================================================

mat gp::calculateKthComponentKernelSubBlock(double amplitude, double length) const {
  mat sub_block = amplitude * exp(time_diff_mat / length);
  sub_block.diag() += kernel_jitter * amplitude;
  return sub_block;
};

void gp::calculateKernelSubBlock() {
  for(uword k = 0; k < K; k++) {
    kernel_sub_block.slice(k) = calculateKthComponentKernelSubBlock(amplitude(k), length(k));
  }
};

// === Parameter updates =======================================================

double gp::muLogDensity(const vec& mu_k, const mat& kernel) const {
  const mat Lower = cholLowerRobust(kernel);
  const vec z = solve(trimatl(Lower), mu_k);
  return -0.5 * ((double) P * std::log(2.0 * M_PI) + 2.0 * accu(log(Lower.diag())) + dot(z, z));
}

void gp::sampleAmplitudeAndLength(uword k) {
  const vec mu_k = mu.col(k);
  const double log_prior_scale = hyper_prior_sd * hyper_prior_sd;
  
  // Target for (log amplitude, log length): N(mu_k; 0, K) times the log-normal 
  // priors, i.e. the log-normal density in the log of the hyperparameter.
  auto log_target = [&](double a, double l, const mat& kernel) {
    return muLogDensity(mu_k, kernel) 
      + pNorm(std::log(a), 0.0, log_prior_scale) 
      + pNorm(std::log(l), 0.0, log_prior_scale);
  };
  
  double current = log_target(amplitude(k), length(k), kernel_sub_block.slice(k));
  
  // Amplitude
  double proposal = amplitude(k) * std::exp(amplitude_proposal_window * randn());
  if(proposal >= hyper_lower && proposal <= hyper_upper) {
    const mat kernel = calculateKthComponentKernelSubBlock(proposal, length(k));
    const double proposed = log_target(proposal, length(k), kernel);
    if(std::log(randu()) < proposed - current) {
      amplitude(k) = proposal;
      kernel_sub_block.slice(k) = kernel;
      current = proposed;
      amplitude_acceptance_count(k)++;
    }
  }
  
  // Length
  proposal = length(k) * std::exp(length_proposal_window * randn());
  if(proposal >= hyper_lower && proposal <= hyper_upper) {
    const mat kernel = calculateKthComponentKernelSubBlock(amplitude(k), proposal);
    const double proposed = log_target(amplitude(k), proposal, kernel);
    if(std::log(randu()) < proposed - current) {
      length(k) = proposal;
      kernel_sub_block.slice(k) = kernel;
      length_acceptance_count(k)++;
    }
  }
}

void gp::sampleNoise(uword k, const mat& component_data) {
  const double n_k = (double) component_data.n_rows;
  const double log_prior_scale = noise_prior_sd * noise_prior_sd;
  const double sum_sq = accu(square(component_data.each_row() - mu.col(k).t()));
  
  auto log_target = [&](double s) {
    return -0.5 * sum_sq / s - 0.5 * n_k * (double) P * std::log(s) 
      + pNorm(std::log(s), 0.0, log_prior_scale);
  };
  
  const double proposal = noise(k) * std::exp(noise_proposal_window * randn());
  if(proposal < hyper_lower || proposal > hyper_upper) {
    return;
  }
  if(std::log(randu()) < log_target(proposal) - log_target(noise(k))) {
    noise(k) = proposal;
    noise_acceptance_count(k)++;
  }
}

void gp::sampleKthComponentParameters(uword k, const umat& members, const uvec& non_outliers) {
  
  const uvec rel_inds = find((members.col(k) == 1) && (non_outliers == 1));
  const uword n_k = rel_inds.n_elem;
  
  if(n_k > 0){
    const mat component_data = X.rows( rel_inds ) ;
    const vec sample_mean = mean(component_data, 0).t();
    
    // Posterior of mu_k: N(K Q^{-1} (n / noise) xbar, K Q^{-1}) with 
    // Q = I + (n / noise) K. K and Q commute, so K Q^{-1} is symmetric.
    const mat& kernel = kernel_sub_block.slice(k);
    const mat Q = I_p + ((double) n_k / noise(k)) * kernel;
    const mat cov_tilde = solve(Q, kernel);
    const vec mu_tilde = ((double) n_k / noise(k)) * (cov_tilde * sample_mean);
    mu.col(k) = rmvnormChol(mu_tilde, cov_tilde);
    
    const bool update_hypers = (samplingCount % sampleHypersFrequency) == 0;
    if(update_hypers) {
      sampleAmplitudeAndLength(k);
      sampleNoise(k, component_data);
    }
  } else {
    // Empty components are drawn from the prior; the kernel must be built from 
    // the new hyperparameters before mu is drawn
    sampleKthComponentHyperParameterPrior(k);
    mu.col(k) = rmvnormChol(zeros<vec>(P), kernel_sub_block.slice(k));
  }
};

void gp::sampleParameters(const arma::umat& members, const arma::uvec& non_outliers) {
  calculateKernelSubBlock();
  for(uword k = 0; k < K; k++) {
    sampleKthComponentParameters(k, members, non_outliers);
  }
  samplingCount++;
  recordHypers();
};

void gp::receiveHyperParametersProposalWindows(vec proposal_windows) {
  if(proposal_windows.n_elem < 3) {
    Rcpp::stop("GP proposal windows must have three entries (amplitude, length, noise).");
  }
  amplitude_proposal_window = proposal_windows[0];
  length_proposal_window = proposal_windows[1];
  noise_proposal_window = proposal_windows[2];
}

// === Likelihood and missing data =============================================

double gp::logLikelihood(arma::uword n, arma::uword k) {
  const arma::uvec& obs_idx = observed_indices(n);
  double ss = 0.0;
  for(uword i = 0; i < obs_idx.n_elem; i++) {
    const uword p = obs_idx(i);
    ss += std::pow(X(n, p) - mu(p, k), 2.0);
  }
  return -0.5 * ss / noise(k) 
    - 0.5 * (double) obs_idx.n_elem * (std::log(2.0 * M_PI) + std::log(noise(k)));
};

arma::vec gp::itemLogLikelihood(arma::uword n) {
  arma::vec ll(K);
  for(uword k = 0; k < K; k++) {
    ll(k) = logLikelihood(n, k);
  }
  return ll;
}

void gp::sampleMissingForObservation(arma::uword n) {
  const uword k = labels(n);
  const arma::uvec& miss_idx = missing_indices(n);
  const double sd = std::sqrt(noise(k));
  for(uword i = 0; i < miss_idx.n_elem; i++) {
    const uword p = miss_idx(i);
    X(n, p) = mu(p, k) + sd * randn();
  }
}

// === Relabelling and predictive checks =======================================

void gp::swapComponents(uword k, uword kprime) {
  mu.swap_cols(k, kprime);
  std::swap(amplitude(k), amplitude(kprime));
  std::swap(length(k), length(kprime));
  std::swap(noise(k), noise(kprime));
  kernel_sub_block.slice(k).swap(kernel_sub_block.slice(kprime));
  std::swap(amplitude_acceptance_count(k), amplitude_acceptance_count(kprime));
  std::swap(length_acceptance_count(k), length_acceptance_count(kprime));
  std::swap(noise_acceptance_count(k), noise_acceptance_count(kprime));
  recordHypers();
}

arma::vec gp::parameters() const {
  return join_cols(vectorise(mu), noise);
}

void gp::setParameters(const arma::vec& theta) {
  if(theta.n_elem != P * K + K) {
    Rcpp::stop("gp: parameter vector has the wrong length.");
  }
  mu = reshape(theta.subvec(0, P * K - 1), P, K);
  noise = theta.subvec(P * K, P * K + K - 1);
}

arma::vec gp::simulate(arma::uword k) const {
  return mu.col(k) + std::sqrt(noise(k)) * randn<arma::vec>(P);
}
