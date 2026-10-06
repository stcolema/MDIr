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

// The complete state of one sampler, enough to copy a particle (see mdi::saveState)
struct mdiState {
  arma::umat labels;
  arma::mat w;
  arma::vec phis, mass;
  double v = 0.0, Z = 0.0;
  std::vector<arma::vec> theta, pooled;
};

class mdi {

public:

  uword N, L, K_max, LC2 = 1, acceptance_count = 0;

  double
    // Normalising constant
    Z = 0.0,

    // Z at the start of the latest sweep (before any phi is updated), the value
    // recorded in the trace of the normalising constant
    Z_start = 0.0,

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
    outlier_types;      // outliers types used

  // For each view, the components that hold no item with an observed label (any
  // index, not necessarily the last ones). Only these are exchanged by the
  // label-swap move; a component with an observed member never moves.
  std::vector<arma::uvec> free_components;

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

  // Rate (without the prior rate) of the conditional of phi_lm
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

  // Draw each phi from its conditional with the strategic latent variable
  // integrated out (see mdiSamplePhiSlice()). Together with the draw of v that
  // follows it at the start of a sweep this is a joint draw of (phi, v) given
  // the weights and labels. It must come before v is redrawn, because the 
  // new phi makes the current v stale.
  void updatePhisSlice();

  // Use updatePhisSlice() (TRUE) or the Gibbs update given v (FALSE)
  bool phi_slice = true;

  // === Tempering =============================================================

  // Inverse temperature. The tempered target is
  //   pi_beta(state) proportional to L(state)^beta P(state),
  // where L is the likelihood of the data given the labels and the component
  // parameters, sum_{l,n} log f_l(x_nl | theta_{l, c_nl}), and P is everything
  // else (the priors, the coupling of the labels through the weights and phis,
  // and the strategic latent variable). beta = 1 is the posterior. Only the 
  // allocation and the component parameters involve L, so only those updates 
  // change with beta; the weights, phis, masses, pooled hyperparameters and the
  // strategic latent variable have the same conditionals at every beta.
  double beta = 1.0;
  void setBeta(double beta_new);

  // log L at the current state: sum over views and items of the log-density of
  // the item at its component. Unchanged by a relabelling that exchanges 
  // component parameters with the labels (updateLabels()).
  double dataLogLikelihood();

  // Copy a state out of, and into, a sampler (used to resample particles). The
  // state holds the labels, weights, phis, masses, strategic latent variable, the
  // component parameters and the pooled hyperparameters; everything else is
  // recomputed. Supported for the densities whose parameters() and
  // setParameters() are complete (G, MVN, C).
  mdiState saveState() const;
  void loadState(const mdiState& state);

  // Replace the state by an exact draw from the prior: masses, phis, weights,
  // component parameters and pooled hyperparameters from their priors, then the
  // labels of every item from p(c | w, phi). The strategic latent variable is
  // drawn from its conditional given the weights and phis. This is a draw from
  // the beta = 0 target.
  void initialiseFromPrior();

  // === Likelihood ============================================================

  // log p(x_n | w, phi, theta) for every item at the current state: the data
  // of item n marginalised over all joint component assignments (and outlier
  // status), which is Z(w * g_n) / Z(w) (see mdiPartition.h). Items with an
  // observed label in a view contribute the density of the data and of that label
  // (the sum is restricted to the observed component), as the label is data.
  arma::vec pointwiseLogLikelihood();

  // The log of g(k, l) = p(x_nl | component k of view l) for every view, with
  // the outlier distribution marginalised where a view has one, in the K_max x L
  // matrix log_g (-Inf beyond K(l)). If use_fixed, a view in which item n has an
  // observed label keeps only that component (and no outlier).
  void componentLogLikelihoods(uword n, bool use_fixed, arma::mat& log_g);

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

  // One full Gibbs sweep. The order is: the phis with the strategic latent
  // variable integrated out (if phi_slice; see updatePhisSlice()), the strategic
  // latent variable, mass, weights and (otherwise) phis given it, then the
  // component parameters
  // (given the labels and the current imputations), then the (label, outlier) 
  // draw followed by the imputation of missing values, and every tenth sweep a
  // label swap move within views.
  void sweep(uword iteration);
  
  // Draw the labels of n items jointly across views from p(c | w, phi) by 
  // enumerating all K_1 x ... x K_L combinations (used for prior predictive draws)
  arma::umat samplePriorLabels(uword n_items) const;
  
  void initialiseMDI();
  void initialiseDatasetL(uword l);

  // === Split-merge ===========================================================

  // Attempts per view and sweep of the sequentially-allocated re-partition of two
  // randomly chosen components (mixtureModel::splitMergeMove); 0 turns the move off. The
  // weights, phis and the other views' labels are held fixed, so Z and the
  // strategic latent variable are unchanged.
  uword split_merge_moves = 0, split_merge_attempts = 0, split_merge_accepts = 0;
  void setSplitMerge(uword moves);
  void updateSplitMerge();
  void updateSplitMergeViewL(uword l);

  // === Joint allocation ======================================================

  // Size of the blocks of views whose labels are redrawn together for each item (0 turns
  // the move off). The labels of one item in the views of a block are drawn from their
  // exact joint conditional given the weights, phis, component parameters and the labels in
  // the other views (mdiSampleJointBlock()); with a block as large as the number of views the
  // item is redrawn across all views at once. It is a block of the Gibbs sweep, so it leaves
  // the target invariant, and it removes the factor (1 + phi) that a move in one view alone pays
  // for breaking agreement with the other views. It is followed by the redraw of any
  // missing values, as the ordinary allocation step is.
  uword joint_block = 0;
  void setJointAllocation(uword block_size);
  void updateJointAllocation();

  // === Label swapping ========================================================

  // Metropolis-Hastings relabelling moves within a view. Swaps the labels,
  // weights and component parameters of two components of that view (see
  // mdiSwapLogRatio() for the acceptance ratio); this improves the alignment of
  // clusters across views.
  void updateLabels();
  void updateLabelsViewL(uword lstar);

};

#endif /* MDI_H */
