// mvn.h
// =============================================================================
// include guard
#ifndef MVN_H
#define MVN_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "density.h"
# include "genericFunctions.h"

using namespace arma ;

// =============================================================================
// mvn class
//
// Multivariate normal components with a conjugate normal-inverse-Wishart prior,
//   Sigma_k ~ IW(scale, nu),  mu_k | Sigma_k ~ N(xi, Sigma_k / kappa).
// xi is the column mean of the observed data and scale = (mean marginal
// variance) / K^{2/P} * I; kappa = 0.01 and nu = P + 2 are fixed (see
// hyperparameterList()).
class mvn : virtual public density
{
public:

  // Parameters and hyperparameters
  double kappa, nu;

  arma::vec xi, cov_log_det;
  arma::mat scale, mu;
  arma::cube cov, cov_inv;

  // Pooling of the scale (see above)
  double scale_shape;
  arma::vec scale_prior_mean;
  
  mvn(arma::uword _K, arma::uvec _labels, arma::mat _X, arma::vec _density_prior = arma::vec());

  virtual ~mvn() { };

  // Data-driven hyperparameters
  arma::vec empiricalMean();
  arma::mat empiricalScaleMatrix();
  void empiricalBayesHyperparameters();

  // Sampling from priors
  void sampleCovPrior();
  void sampleMuPrior();
  void sampleFromPriors() override;
  void updatePooledHyperparameters(const arma::uvec& occupied) override;
  arma::vec pooledHyperparameters() const override;
  void setPooledHyperparameters(const arma::vec& pooled) override;

  void sampleKthComponentParameters(uword k, const umat& members, const uvec& non_outliers) override;

  // Cache the inverse and log determinant of every component covariance
  void matrixCombinations();

  // Missing values
  void sampleMissingForObservation(arma::uword n) override;

  // Likelihood of the observed entries
  arma::vec itemLogLikelihood(arma::uword n) override;
  double logLikelihood(arma::uword n, arma::uword k) override;

  void swapComponents(uword k, uword kprime) override;

  // Layout: mu (P x K), then cov (P x P x K), each column-major
  arma::vec parameters() const override;
  void setParameters(const arma::vec& theta) override;
  arma::vec simulate(arma::uword k) const override;

  Rcpp::List hyperparameterList() const override;
};

#endif /* MVN_H */
