// mdi.h
// =============================================================================
// include guard
#ifndef MDI_H
#define MDI_H

// =============================================================================
// included dependencies
# include "logLikelihoods.h"
# include "genericFunctions.h"
# include "mdiPartition.h"
# include "mixtureModel.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// MDI class
//
// The joint model for the L views is
//
//   p(c | w, phi) = prod_n [ prod_l w(c_nl, l) prod_{l<m} (1 + phi_lm 1[c_nl = c_nm]) ] / Z^N
//
// with w(k, l) ~ Gamma(mass_l / K_l, w_rate_prior), phi_lm ~ Gamma(phi_shape_prior,
// phi_rate_prior), mass_l ~ Gamma(mass_shape_prior, mass_rate_prior). The
// intractable Z^{-N} is handled with the strategic latent variable v ~
// Gamma(N, Z) (Kirk et al., 2012), so every parameter has a tractable full
// conditional.

class mdi {

public:

  uword N, L, K_max, LC2 = 1, acceptance_count = 0;

  double
    // Normalising constant
    Z = 0.0,

    // Strategic latent variable
    v = 0.0,

    // Prior hyperparameters for view mass components (log-scale random walk
    // is used to update them)
    mass_proposal_sd = 0.3,
    mass_shape_prior = 2.0,
    mass_rate_prior = 0.1,

    // Prior hyperparameters for component weights
    w_rate_prior = 2.0,

    // Prior hyperparameters for MDI phi parameters
    phi_shape_prior = 2.0,
    phi_rate_prior = 0.2,

    // Model fit
    complete_likelihood = 0.0,
    observed_likelihood = 0.0;

  // Metropolis-Hastings acceptance tallies for the mass parameters
  arma::vec mass_acceptance_count;

  arma::uvec
    K,                  // Number of clusters in each dataset
    mixture_types,      // mixture types used
    outlier_types,      // outliers types used
    K_unfixed,          // Number of components not fixed
    K_fixed;            // Number of components fixed (i.e. at least one member has an observed label)

  // Options of the density-level priors (see resolveDensityPrior())
  arma::vec density_prior;
  
  arma::vec phis,
    mass,
    complete_likelihood_vec,
    observed_likelihood_vec;

  arma::umat
    labels,

    // Map from a pair of views to the index of their phi
    phi_map,

    // Class membership in each dataset
    N_k,

    // Indicator matrix for item n being an outlier in dataset l
    outliers,

    // Indicator matrix for item n being well-described by its component
    // in dataset l
    non_outliers,

    fixed;

  // The weights in each dataset
  arma::mat w;

  // Connected sums of the partition-sum recursion and the phis they were built
  // from (see mdiPartition.h)
  std::vector<double> partition_tables;
  arma::vec partition_tables_phis;

  // Cube of cluster members
  arma::ucube members;

  // The data can have varying numbers of columns
  arma::field<arma::mat> X;

  // The collection of mixtures
  std::vector< std::unique_ptr<mixtureModel> > mixtures;

  mdi(
    arma::field<arma::mat> _X,
    uvec _mixture_types,
    uvec _outlier_types,
    arma::uvec _K,
    arma::umat _labels,
    arma::umat _fixed,
    arma::vec _prior = arma::vec(),
    arma::vec _density_prior = arma::vec()
  ) ;

  virtual ~mdi() { };
  
  // Set the MDI-level prior hyperparameters from a vector ordered as
  // (mass_shape, mass_rate, weight_rate, phi_shape, phi_rate); an empty vector 
  // keeps the defaults.
  void setPrior(const arma::vec& prior);

  // === Normalising constant, weights and phis ================================

  // The L x L symmetric matrix of phis
  arma::mat phiMatrix() const;

  // Refresh the cached connected sums if the phis have changed since they were
  // computed
  void refreshPartitionTables();

  void updateNormalisingConstant();
  void sampleStrategicLatentVariable();

  // Rate (without the prior rate) of the conditional of w(k, l), and of phi_lm
  double calcWeightRate(uword lstar, uword kstar) const;
  double calcPhiRate(uword l, uword m) const;

  void updateWeights();
  void updateWeightsViewL(uword lstar);

  void updateMassParameters();
  void updateMassParameterViewL(uword lstar);

  // Shape of the Gamma(phi_shape_prior + r, phi_rate_prior + rate) mixture
  // over the phi posterior. Returns the log-weights over r = 0, ..., N_lm.
  arma::vec calculatePhiShapeMixtureWeights(uword N_lm, double rate) const;
  uword samplePhiShape(uword N_lm, double rate) const;
  void updatePhis();

  // === Allocations ===========================================================

  // Recompute the membership indicators and counts of view l from its labels
  void refreshMembersViewL(uword l);
  void updateAllocation();
  void updateAllocationViewL(uword l);

  // === Initialisation ========================================================

  void initialiseMixtures();
  void sampleFromPriors();
  void sampleFromLocalPriors();
  void sampleFromGlobalPriors();
  vec samplePhiPrior(uword n_phis);
  double sampleWeightPrior(uword l);
  vec sampleMassPrior();

  // log(1 + phi) for the items that share a label with item n in another view
  mat calculateUpweights(uword l) const;

  // One full Gibbs sweep. The order is: strategic latent variable, mass, 
  // weights and phis (given the current labels), then the component parameters
  // (given the labels and the current imputations), then the (label, outlier) 
  // draw followed by the imputation of missing values, and every tenth sweep a
  // label swap move within views.
  void sweep(uword iteration);
  
  // Draw the labels of n items jointly across views from p(c | w, phi) by 
  // enumerating all K_1 x ... x K_L combinations (used for prior predictive draws)
  arma::umat samplePriorLabels(uword n_items) const;
  
  void initialiseMDI();
  void initialiseDatasetL(uword l);

  // === Label swapping ========================================================

  // Metropolis-Hastings relabelling moves within a view. Swaps the labels,
  // weights and component parameters of two components; this improves the
  // alignment of clusters across views.
  double calcScore(uword lstar, const arma::umat& c) const;
  void updateLabels();
  void updateLabelsViewL(uword lstar);

};

#endif /* MDI_H */
