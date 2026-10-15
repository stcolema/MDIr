// gp.h
// =============================================================================
// include guard
#ifndef GP_H
#define GP_H

// =============================================================================
// included dependencies
# include "density.h"

using namespace arma ;

// Parameters (alpha, beta) of the inverse-gamma distribution with
// P(lambda < lower) = tail and P(lambda > upper) = tail. Used to keep the length 
// scale of a Gaussian process within what the measurement grid can resolve.
arma::vec calibrateInverseGamma(double lower, double upper, double tail = 0.01);

// =============================================================================
// gp class
//
// Each component has a mean function over the P equally spaced measurements,
//   mu_k = xi + f_k,   f_k ~ N(0, K_k),   K_k[i, j] = a_k exp(-(i - j)^2 / (2 lambda_k^2)),
// with xi the column means of the observed data, and items are noisy draws 
// x_n ~ N(mu_k, s2_k I). The spacing between adjacent measurements is one unit.
//
// Priors (following the advice not to leave Gaussian process hyperparameters
// free to fit noise, in particular length scales shorter than the grid):
//  * lambda_k has a hard floor at gp_min_length (default: the grid spacing, 1)
//    and an inverse-gamma prior calibrated so that 1% of its mass lies below
//    the floor and 1% above the extent of the grid, P - 1 (Betancourt, 2017, 
//    "Robust Gaussian Process Modeling", Stan case study). A length scale below
//    the spacing makes K_k diagonal, so f_k is white noise that cannot be 
//    told apart from the measurement noise and the mean function fits noise.
//  * log a_k ~ N(m_a, s_a^2) and log s2_k ~ N(m_n, s_n^2), a population shared
//    by the components (partial pooling), with m ~ N(log v, gp_center_sd^2) 
//    centred on the average data variance v and s ~ half-normal(gp_pool_sd_scale)
//    (Gelman, 2006). With gp_pool = 0 the population is fixed at (log v, 1).
//
// Updates: (a_k, lambda_k) by random-walk Metropolis on the log scale with
// target N(f_k; 0, K_k) times the priors, s2_k by Metropolis given f_k and the
// data, plus a non-centred amplitude move that holds f / sqrt(a) fixed; (m, s) given the occupied components' log a_k / log s2_k, marginalising
// m in the update of s to avoid the funnel.
class gp : virtual public density
{
public:
  
  uword
    // Component hyperparameters are updated every `sampleHypersFrequency` iterations
    sampleHypersFrequency = 5, 
    samplingCount = 0;
  
  bool pool = true;
  
  double
    // Data-scale centre for the amplitude and noise populations (log scale)
    log_variance_centre = 0.0,
    center_sd = 2.0,
    pool_sd_scale = 1.0,
    
    // Length scale range and inverse-gamma prior
    min_length = 1.0,
    max_length = 10.0,
    length_shape = 1.0,
    length_rate = 1.0,
    
    // Numerical limits on amplitude and noise
    hyper_lower = 1e-10,
    hyper_upper = 1e10,
    kernel_jitter = 1e-8,
    
    // Populations of log amplitude and log noise: mean and sd
    m_amp = 0.0, s_amp = 1.0, m_noise = 0.0, s_noise = 1.0,
    
    // Random-walk standard deviations on the log scale
    amplitude_proposal_window = 0.25,
    length_proposal_window = 0.25,
    noise_proposal_window = 0.15,
    pool_proposal_window = 0.5;
  
  uvec noise_acceptance_count,
    length_acceptance_count,
    amplitude_acceptance_count,
    noise_attempt_count,
    length_attempt_count,
    amplitude_attempt_count;
  
  vec amplitude, length, noise, xi;
  mat mu, I_p, time_diff_mat;
  cube kernel_sub_block;
  
  gp(arma::uword _K, arma::uvec _labels, arma::mat _X, arma::vec _density_prior = arma::vec());
  
  virtual ~gp() { };
  
  // Priors
  double sampleTruncatedLogNormal(double m, double s) const;
  double sampleLengthPrior() const;
  void sampleKthComponentHyperParameterPrior(uword k);
  void sampleFromPriors() override;
  
  // Log-prior densities in the log of the hyperparameter (Jacobians included)
  double logPriorLogAmplitude(double log_a) const;
  double logPriorLogNoise(double log_s2) const;
  double logPriorLogLength(double length) const;
  
  // Kernel
  mat calculateKthComponentKernelSubBlock(double amplitude, double length) const;
  void calculateKernelSubBlock();
  
  void sampleKthComponentParameters(uword k, const umat& members, const uvec& non_outliers) override;
  void sampleParameters(const arma::umat& members, const arma::uvec& non_outliers) override;
  
  void updatePooledHyperparameters(const arma::uvec& occupied) override;
  arma::vec pooledHyperparameters() const override;
  void updatePopulation(const arma::vec& y, double& m, double& s) const;
  
  void receiveHyperParametersProposalWindows(vec proposal_windows) override;
  
  // Log-density of f_k under N(0, kernel)
  double muLogDensity(const vec& f_k, const mat& kernel) const;
  void sampleAmplitudeAndLength(uword k);
  void sampleAmplitudeNonCentred(uword k, const mat& component_data);
  void sampleNoise(uword k, const mat& component_data);
  
  void sampleMissingForObservation(arma::uword n) override;
  
  arma::vec itemLogLikelihood(arma::uword n) override;
  double logLikelihood(arma::uword n, arma::uword k) override;
  
  void swapComponents(uword k, uword kprime) override;
  
  // Layout: mu (P x K), then noise (K)
  arma::vec parameters() const override;
  void setParameters(const arma::vec& theta) override;
  arma::vec simulate(arma::uword k) const override;
  
  Rcpp::List hyperparameterList() const override;
  
private:
  void recordHypers();
};

#endif /* GP_H */
