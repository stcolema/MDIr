// gaussian.h
// =============================================================================
// include guard
#ifndef GAUSSIAN_H
#define GAUSSIAN_H

// =============================================================================
// included dependencies
# include "density.h"

using namespace arma ;

// =============================================================================
// gaussian class
//
// Independent normal features within each component (diagonal covariance) with
// a conjugate normal-inverse-gamma prior per feature,
//   sigma2_kp ~ InvGamma(nu / 2, scale_p / 2),  mu_kp | sigma2_kp ~ N(xi_p, sigma2_kp / kappa).
class gaussian : virtual public density
{
public:

  // Parameters and hyperparameters
  double kappa = 0.01, nu = 3.0;

  arma::vec xi, scale;

  // Component parameters, P x K. Note these are variances, not standard
  // deviations.
  arma::mat mu, variances, precisions, log_precisions;

  gaussian(arma::uword _K, arma::uvec _labels, arma::mat _X);

  virtual ~gaussian() { };

  // Data-driven hyperparameters
  arma::vec empiricalMean();
  arma::vec empiricalScaleVector();
  void empiricalBayesHyperparameters();

  // Sampling from priors
  void sampleVariancePrior();
  void sampleMuPrior();
  void sampleFromPriors() override;

  void sampleKthComponentParameters(uword k, const umat& members, const uvec& non_outliers) override;

  void sampleMissingForObservation(arma::uword n) override;

  arma::vec itemLogLikelihood(arma::uword n) override;
  double logLikelihood(arma::uword n, arma::uword k) override;

  void swapComponents(uword k, uword kprime) override;

  // Layout: mu (P x K), then variances (P x K), each column-major
  arma::vec parameters() const override;
  void setParameters(const arma::vec& theta) override;
  arma::vec simulate(arma::uword k) const override;

  Rcpp::List hyperparameterList() const override;

private:
  void setKthVariance(uword k, uword p, double variance);
};

#endif /* GAUSSIAN_H */
