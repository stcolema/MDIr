// gp.cpp
// =============================================================================
# include "logLikelihoods.h"
# include "gp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// Inverse-gamma(alpha, beta) with P(lambda < lower) = P(lambda > upper) = tail.
// If lambda ~ InvGamma(alpha, beta) then 1 / lambda ~ Gamma(alpha, rate = beta), so
//   P(lambda < lower) = P(Gamma > 1 / lower),  P(lambda > upper) = P(Gamma < 1 / upper).
// For a given alpha, beta is fixed by the lower tail (P(lambda < lower) decreases
// as beta increases, since lambda tends to be larger); the upper tail then
// decreases in alpha, which is bisected.
arma::vec calibrateInverseGamma(double lower, double upper, double tail) {
  if(!(lower > 0.0) || !(upper > lower) || !(tail > 0.0 && tail < 0.5)) {
    Rcpp::stop("calibrateInverseGamma: need 0 < lower < upper and 0 < tail < 0.5.");
  }
  auto beta_for = [&](double alpha) {
    double lo = 1e-8, hi = 1e8;
    for(int it = 0; it < 200; it++) {
      const double mid = std::sqrt(lo * hi);
      // P(lambda < lower) = P(Gamma(alpha, rate = mid) > 1 / lower)
      const double p = R::pgamma(1.0 / lower, alpha, 1.0 / mid, 0, 0);
      if(p > tail) lo = mid; else hi = mid;
    }
    return std::sqrt(lo * hi);
  };
  auto upper_tail = [&](double alpha, double beta) {
    return R::pgamma(1.0 / upper, alpha, 1.0 / beta, 1, 0);
  };
  double a_lo = 0.05, a_hi = 500.0;
  for(int it = 0; it < 200; it++) {
    const double mid = std::sqrt(a_lo * a_hi);
    if(upper_tail(mid, beta_for(mid)) > tail) a_lo = mid; else a_hi = mid;
  }
  const double alpha = std::sqrt(a_lo * a_hi);
  return arma::vec({alpha, beta_for(alpha)});
}

gp::gp(arma::uword _K, arma::uvec _labels, arma::mat _X, arma::vec _density_prior) : 
  density(_K, _labels, _X, _density_prior) 
{
  if(P < 3) {
    Rcpp::stop("Gaussian process views need at least three measurements per item.");
  }
  
  pool = density_prior(1) > 0.5;
  min_length = density_prior(2);
  center_sd = density_prior(3);
  pool_sd_scale = density_prior(4);
  
  amplitude.ones(K);
  length.ones(K);
  noise.ones(K);
  mu.zeros(P, K);
  kernel_sub_block.zeros(P, P, K);
  I_p = eye(P, P);
  
  // -(i - j)^2 / 2, so that K = amplitude * exp(time_diff_mat / length^2)
  time_diff_mat.zeros(P, P);
  for(uword ii = 0; ii < P; ii++) {
    for(uword jj = ii + 1; jj < P; jj++) {
      time_diff_mat(ii, jj) = - 0.5 * std::pow((double) (jj - ii), 2.0);
      time_diff_mat(jj, ii) = time_diff_mat(ii, jj);
    }
  }
  
  // Length scale: floor at min_length and inverse-gamma prior with 1% of mass 
  // below the floor and 1% above the extent of the grid (at least twice the floor)
  const double extent = std::max((double) (P - 1), 2.0 * min_length);
  const arma::vec ig = calibrateInverseGamma(min_length, extent, 0.01);
  length_shape = ig(0);
  length_rate = ig(1);
  max_length = 10.0 * extent;
  
  // Data-driven location: column means and the average data variance
  xi = sampleMeanRobust(X);
  const arma::mat global_cov = computeCovarianceRobust(X);
  log_variance_centre = std::log(std::max(arma::accu(global_cov.diag()) / (double) P, 1e-12));
  m_amp = m_noise = log_variance_centre;
  s_amp = s_noise = 1.0;
  
  noise_acceptance_count.zeros(K);
  length_acceptance_count.zeros(K);
  amplitude_acceptance_count.zeros(K);
  noise_attempt_count.zeros(K);
  length_attempt_count.zeros(K);
  amplitude_attempt_count.zeros(K);
  
  n_param = P + 3;
  
  hypers.zeros(3 * K);
  acceptance_count.zeros(3 * K);
  acceptance_attempts.zeros(3 * K);
  
  identifyMissingValues();
  initializeMissingValues();
};

Rcpp::List gp::hyperparameterList() const {
  return Rcpp::List::create(
    Rcpp::Named("xi") = xi,
    Rcpp::Named("log_variance_centre") = log_variance_centre,
    Rcpp::Named("pooling") = pool,
    Rcpp::Named("center_sd") = center_sd,
    Rcpp::Named("pool_sd_scale") = pool_sd_scale,
    Rcpp::Named("min_length") = min_length,
    Rcpp::Named("max_length") = max_length,
    Rcpp::Named("length_prior_shape") = length_shape,
    Rcpp::Named("length_prior_rate") = length_rate,
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
  acceptance_attempts.subvec(0, K - 1) = amplitude_attempt_count;
  acceptance_attempts.subvec(K, 2 * K - 1) = length_attempt_count;
  acceptance_attempts.subvec(2 * K, 3 * K - 1) = noise_attempt_count;
}

// === Priors ==================================================================

double gp::sampleTruncatedLogNormal(double m, double s) const {
  double x = 0.0;
  for(int attempt = 0; attempt < 1000; attempt++) {
    x = std::exp(m + s * randn());
    if(x >= hyper_lower && x <= hyper_upper) {
      return x;
    }
  }
  return std::min(std::max(x, hyper_lower), hyper_upper);
}

// Inverse-gamma prior on the length scale restricted to [min_length, max_length]
double gp::sampleLengthPrior() const {
  double x = 0.0;
  for(int attempt = 0; attempt < 1000; attempt++) {
    x = 1.0 / rGamma(length_shape, length_rate);
    if(x >= min_length && x <= max_length) {
      return x;
    }
  }
  return std::min(std::max(x, min_length), max_length);
}

double gp::logPriorLogAmplitude(double log_a) const {
  return -0.5 * std::pow((log_a - m_amp) / s_amp, 2.0);
}

double gp::logPriorLogNoise(double log_s2) const {
  return -0.5 * std::pow((log_s2 - m_noise) / s_noise, 2.0);
}

// Density of log(lambda) when lambda ~ InvGamma(shape, rate): the inverse-gamma
// log-density plus log(lambda) for the change of variables
double gp::logPriorLogLength(double lambda) const {
  return -length_shape * std::log(lambda) - length_rate / lambda;
}

void gp::sampleKthComponentHyperParameterPrior(uword k) {
  amplitude(k) = sampleTruncatedLogNormal(m_amp, s_amp);
  noise(k) = sampleTruncatedLogNormal(m_noise, s_noise);
  length(k) = sampleLengthPrior();
  kernel_sub_block.slice(k) = calculateKthComponentKernelSubBlock(amplitude(k), length(k));
};

void gp::sampleFromPriors() {
  if(pool) {
    m_amp = log_variance_centre + center_sd * randn();
    m_noise = log_variance_centre + center_sd * randn();
    s_amp = pool_sd_scale * std::abs(randn());
    s_noise = pool_sd_scale * std::abs(randn());
    s_amp = std::max(s_amp, 1e-6);
    s_noise = std::max(s_noise, 1e-6);
  } else {
    m_amp = m_noise = log_variance_centre;
    s_amp = s_noise = 1.0;
  }
  for(uword k = 0; k < K; k++) {
    sampleKthComponentHyperParameterPrior(k);
    mu.col(k) = xi + rmvnormChol(zeros<vec>(P), kernel_sub_block.slice(k));
  }
  recordHypers();
};

// === Kernel ==================================================================

mat gp::calculateKthComponentKernelSubBlock(double amplitude, double length) const {
  mat sub_block = amplitude * exp(time_diff_mat / (length * length));
  sub_block.diag() += kernel_jitter * amplitude;
  return sub_block;
};

void gp::calculateKernelSubBlock() {
  for(uword k = 0; k < K; k++) {
    kernel_sub_block.slice(k) = calculateKthComponentKernelSubBlock(amplitude(k), length(k));
  }
};

// === Parameter updates =======================================================

double gp::muLogDensity(const vec& f_k, const mat& kernel) const {
  const mat Lower = cholLowerRobust(kernel);
  const vec z = solve(trimatl(Lower), f_k);
  return -0.5 * ((double) P * std::log(2.0 * M_PI) + 2.0 * accu(log(Lower.diag())) + dot(z, z));
}

// A random-walk scale drawn from {0.4, 1, 2.5} times the base window. A mixture
// of symmetric random walks is itself symmetric, and it copes with amplitudes 
// that are either tightly or loosely identified by the data.
static double mixedWindow(double window) {
  const double u = randu();
  return window * (u < 1.0 / 3.0 ? 0.4 : (u < 2.0 / 3.0 ? 1.0 : 2.5));
}

void gp::sampleAmplitudeAndLength(uword k) {
  const vec f_k = mu.col(k) - xi;
  
  // Target for (log amplitude, log length): N(f_k; 0, K) times the priors, each
  // expressed as a density in the log of the hyperparameter.
  auto log_target = [&](double a, double l, const mat& kernel) {
    return muLogDensity(f_k, kernel) + logPriorLogAmplitude(std::log(a)) + logPriorLogLength(l);
  };
  
  double current = log_target(amplitude(k), length(k), kernel_sub_block.slice(k));
  
  // Amplitude
  double proposal = amplitude(k) * std::exp(mixedWindow(amplitude_proposal_window) * randn());
  amplitude_attempt_count(k)++;
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
  
  // Length; proposals below the floor (or above the ceiling) have prior zero
  proposal = length(k) * std::exp(length_proposal_window * randn());
  length_attempt_count(k)++;
  if(proposal >= min_length && proposal <= max_length) {
    const mat kernel = calculateKthComponentKernelSubBlock(amplitude(k), proposal);
    const double proposed = log_target(amplitude(k), proposal, kernel);
    if(std::log(randu()) < proposed - current) {
      length(k) = proposal;
      kernel_sub_block.slice(k) = kernel;
      length_acceptance_count(k)++;
    }
  }
}

// Non-centred amplitude move. With e = f / sqrt(a) the prior of e does not depend 
// on the amplitude, so proposing a' while holding e fixed (f' = sqrt(a' / a) f) is a
// Metropolis step on the reparametrised target whose ratio involves only the
// likelihood of the component data and the prior of a. It moves the amplitude
// when the data barely constrain the mean function, which is where the centred
// update (given f) mixes slowly; alternating the two is an interweaving strategy
// (Yu and Meng, 2011).
void gp::sampleAmplitudeNonCentred(uword k, const mat& component_data) {
  const double n_k = (double) component_data.n_rows;
  const vec x_bar = mean(component_data, 0).t();
  const double proposal = amplitude(k) * std::exp(mixedWindow(amplitude_proposal_window) * randn());
  amplitude_attempt_count(k)++;
  if(proposal < hyper_lower || proposal > hyper_upper) {
    return;
  }
  const double scale = std::sqrt(proposal / amplitude(k));
  const vec f_current = mu.col(k) - xi;
  const vec mu_proposed = xi + scale * f_current;
  
  // The data enter through n ||xbar - mu||^2 / (2 noise)
  const double log_ratio = 
    -0.5 * n_k * (accu(square(x_bar - mu_proposed)) - accu(square(x_bar - mu.col(k)))) / noise(k)
    + logPriorLogAmplitude(std::log(proposal)) - logPriorLogAmplitude(std::log(amplitude(k)));
  
  if(std::log(randu()) < log_ratio) {
    amplitude(k) = proposal;
    mu.col(k) = mu_proposed;
    kernel_sub_block.slice(k) = calculateKthComponentKernelSubBlock(amplitude(k), length(k));
    amplitude_acceptance_count(k)++;
  }
}

void gp::sampleNoise(uword k, const mat& component_data) {
  const double n_k = (double) component_data.n_rows;
  const double sum_sq = accu(square(component_data.each_row() - mu.col(k).t()));
  
  auto log_target = [&](double s) {
    return -0.5 * sum_sq / s - 0.5 * n_k * (double) P * std::log(s) + logPriorLogNoise(std::log(s));
  };
  
  const double proposal = noise(k) * std::exp(noise_proposal_window * randn());
  noise_attempt_count(k)++;
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
    const vec sample_mean = mean(component_data, 0).t() - xi;
    
    // Posterior of f_k: N(K Q^{-1} (n / noise) xbar, K Q^{-1}) with 
    // Q = I + (n / noise) K. K and Q commute, so K Q^{-1} is symmetric.
    const mat& kernel = kernel_sub_block.slice(k);
    const mat Q = I_p + ((double) n_k / noise(k)) * kernel;
    const mat cov_tilde = solve(Q, kernel);
    const vec mu_tilde = ((double) n_k / noise(k)) * (cov_tilde * sample_mean);
    mu.col(k) = xi + rmvnormChol(mu_tilde, cov_tilde);
    
    if((samplingCount % sampleHypersFrequency) == 0) {
      sampleAmplitudeAndLength(k);
      sampleAmplitudeNonCentred(k, component_data);
      sampleNoise(k, component_data);
    }
  } else {
    // Empty components are drawn from the prior; the kernel must be built from 
    // the new hyperparameters before mu is drawn
    sampleKthComponentHyperParameterPrior(k);
    mu.col(k) = xi + rmvnormChol(zeros<vec>(P), kernel_sub_block.slice(k));
  }
};

void gp::sampleParameters(const arma::umat& members, const arma::uvec& non_outliers) {
  calculateKernelSubBlock();
  density::sampleParameters(members, non_outliers);
  samplingCount++;
  recordHypers();
};

// Update the mean m and sd s of a population of log hyperparameters y, with
// m ~ N(log_variance_centre, center_sd^2) and s ~ half-normal(pool_sd_scale). 
// s is updated by Metropolis on log s using the density of y with m integrated 
// out; m is then drawn exactly given s.
void gp::updatePopulation(const arma::vec& y, double& m, double& s) const {
  const double m0 = log_variance_centre, t2 = center_sd * center_sd;
  const double n = (double) y.n_elem;
  
  if(y.n_elem == 0) {
    m = m0 + center_sd * randn();
    s = std::max(pool_sd_scale * std::abs(randn()), 1e-6);
    return;
  }
  
  const double y_bar = arma::mean(y);
  const double ss_within = arma::accu(arma::square(y - y_bar));
  
  auto log_target = [&](double log_s) {
    const double s2 = std::exp(2.0 * log_s);
    // density of y given s (m integrated out), half-normal prior, Jacobian of log s
    return -(n - 1.0) * log_s - 0.5 * std::log(s2 + n * t2)
      - 0.5 * (ss_within / s2 + n * std::pow(y_bar - m0, 2.0) / (s2 + n * t2))
      - 0.5 * s2 / (pool_sd_scale * pool_sd_scale) + log_s;
  };
  
  const double current = std::log(s);
  const double proposal = current + pool_proposal_window * randn();
  if(std::log(randu()) < log_target(proposal) - log_target(current)) {
    s = std::exp(proposal);
  }
  s = std::max(s, 1e-6);
  
  const double precision = 1.0 / t2 + n / (s * s);
  const double post_mean = (m0 / t2 + arma::accu(y) / (s * s)) / precision;
  m = post_mean + randn() / std::sqrt(precision);
}

void gp::updatePooledHyperparameters(const arma::uvec& occupied) {
  if(!pool) {
    return;
  }
  arma::vec log_amp(occupied.n_elem), log_noise(occupied.n_elem);
  for(uword i = 0; i < occupied.n_elem; i++) {
    log_amp(i) = std::log(amplitude(occupied(i)));
    log_noise(i) = std::log(noise(occupied(i)));
  }
  updatePopulation(log_amp, m_amp, s_amp);
  updatePopulation(log_noise, m_noise, s_noise);
}

arma::vec gp::pooledHyperparameters() const {
  return arma::vec({m_amp, s_amp, m_noise, s_noise});
}

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
  std::swap(amplitude_attempt_count(k), amplitude_attempt_count(kprime));
  std::swap(length_attempt_count(k), length_attempt_count(kprime));
  std::swap(noise_attempt_count(k), noise_attempt_count(kprime));
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
