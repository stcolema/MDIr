// gp.h
// =============================================================================
// include guard
#ifndef GP_H
#define GP_H

// =============================================================================
// included dependencies
# include "density.h"

using namespace arma ;

// =============================================================================
// gp class
//
// Each component has a mean function mu_k over the P (equally spaced) features
// with a Gaussian process prior, mu_k ~ N(0, K(amplitude_k, length_k)), where 
// K is the squared exponential kernel
//   K_ij = amplitude * exp(-(i - j)^2 / (2 * length)),
// and items are noisy draws x_n ~ N(mu_k, noise_k I). The hyperparameters
// (amplitude, length, noise) have independent log-normal priors, log(.) ~ N(0, 1),
// and are updated with random-walk Metropolis steps on the log scale. The 
// amplitude and length are updated conditional on mu_k (target 
// N(mu_k; 0, K(amplitude, length))) and the noise conditional on mu_k and the 
// component data.
class gp : virtual public density
{
public:
  
  uword
    // Hyperparameters are updated every `sampleHypersFrequency` iterations
    sampleHypersFrequency = 5, 
    samplingCount = 0;
  
  double
    // Log-normal prior standard deviations on the log scale
    hyper_prior_sd = 1.0,
    noise_prior_sd = 1.0,
    
    // Hyperparameters are restricted to [lower, upper] to keep the kernel
    // numerically well-behaved
    hyper_lower = 1e-6,
    hyper_upper = 1e6,
    
    // Relative jitter added to the kernel diagonal
    kernel_jitter = 1e-8,
    
    // Random-walk standard deviations on the log scale
    amplitude_proposal_window = 0.25,
    length_proposal_window = 0.25,
    noise_proposal_window = 0.15;
  
  uvec noise_acceptance_count,
    length_acceptance_count,
    amplitude_acceptance_count;
  
  vec amplitude, length, noise;
  mat mu, I_p, time_diff_mat;
  cube kernel_sub_block;
  
  gp(arma::uword _K, arma::uvec _labels, arma::mat _X);
  
  virtual ~gp() { };
  
  // Priors
  double sampleHyperPrior(double sd) const;
  void sampleKthComponentHyperParameterPrior(uword k);
  void sampleFromPriors() override;
  
  // Kernel
  mat calculateKthComponentKernelSubBlock(double amplitude, double length) const;
  void calculateKernelSubBlock();
  
  void sampleKthComponentParameters(uword k, const umat& members, const uvec& non_outliers) override;
  void sampleParameters(const arma::umat& members, const arma::uvec& non_outliers) override;
  
  void receiveHyperParametersProposalWindows(vec proposal_windows) override;
  
  // Log-density of mu_k under N(0, kernel)
  double muLogDensity(const vec& mu_k, const mat& kernel) const;
  void sampleAmplitudeAndLength(uword k);
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
