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
// The scale is pooled across components: scale_p ~ Gamma(a, a / c_p) with shape 
// a = density_prior[0] and prior mean c_p = (mean marginal variance) / K^{2/P}
// (a = 0 fixes scale at c). Then 
//   scale_p | sigma2_occupied ~ Gamma(a + n_occ nu / 2, a / c_p + sum_k 1 / (2 sigma2_kp)).
class gaussian : virtual public density
{
public:

  // Parameters and hyperparameters
  double kappa = 0.01, nu = 3.0;

  // Log-uniform hyperprior on kappa on [kappa_min, kappa_max], shared by the
  // components and features; sampled only if kappa_sampled. Given the occupied
  // components, kappa | . ~ Gamma(P n_occ / 2, sum_kp (mu_kp - xi_p)^2 / (2 sigma2_kp))
  // truncated to the range.
  double kappa_min = 0.0, kappa_max = 0.0;
  bool kappa_sampled = false;
  void updateKappa(const arma::uvec& occupied);

  arma::vec xi, scale;

  // Component parameters, P x K. Note these are variances, not standard
  // deviations.
  arma::mat mu, variances, precisions, log_precisions;

  double scale_shape;
  arma::vec scale_prior_mean;
  
  gaussian(arma::uword _K, arma::uvec _labels, arma::mat _X, arma::vec _density_prior = arma::vec());

  virtual ~gaussian() { };

  // Data-driven hyperparameters
  arma::vec empiricalMean();
  arma::vec empiricalScaleVector();
  void empiricalBayesHyperparameters();

  // Sampling from priors
  void sampleVariancePrior();
  void sampleMuPrior();
  void sampleFromPriors() override;
  void updatePooledHyperparameters(const arma::uvec& occupied) override;
  arma::vec pooledHyperparameters() const override;
  void setPooledHyperparameters(const arma::vec& pooled) override;

  void sampleKthComponentParameters(uword k, const umat& members, const uvec& non_outliers) override;

  void sampleMissingForObservation(arma::uword n) override;

  arma::vec itemLogLikelihood(arma::uword n) override;
  double logLikelihood(arma::uword n, arma::uword k) override;

  // Collapsed normal-inverse-gamma marginal likelihood (independent measurements)
  bool hasCollapsedMarginal() const override { return true; }
  collapsedStats emptyStats() const override;
  void addItemToStats(collapsedStats& st, arma::uword n) const override;
  double logMarginalLikelihood(const collapsedStats& st, double beta) const override;

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
