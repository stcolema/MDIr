// mdi.h
// =============================================================================
// include guard
#ifndef MDI_H
#define MDI_H

// =============================================================================
// included dependencies
# include "logLikelihoods.h"
# include "genericFunctions.h"
# include "mixtureModel.h"

using namespace arma ;

// =============================================================================
// MDI class

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;


class mdiModelAlt {
  
private:
  
public:
  
  bool use_log_norm_proposal = true;
  
  uword N, L, K_max, K_prod, K_to_the_L, n_combinations, LC2 = 1;

  double 
    // Normalising constant
    Z = 0.0,
    
    // Strategic latent variable
    v = 0.0,
    
    // Prior hyperparameters for view mass components
    mass_proposal_window = 0.025,
    mass_shape_prior = 2.0,
    mass_rate_prior = 0.1,
  
    // Prior hyperparameters for component weights
    // w_shape_prior = 2.0,
    w_rate_prior = 2.0,
    
    // Prior hyperparameters for MDI phi parameters
    phi_shape_prior = 1.0,
    phi_rate_prior = 0.1,
    
    // Model fit
    complete_likelihood = 0.0;
  
  
  arma::uvec 
    K,                  // Number of clusters in each dataset
    one_to_K,           // [0, 1, ..., K]
    one_to_L,           // [0, 1, ..., L] 
    KL_powers,          // K to the power of the members of one_to_L
    rows_to_shed,       // rows initialised assuming symmetric K that are shed
    mixture_types,      // mixture types used
    outlier_types,      // outliers types used 
    K_unfixed,          // Number of components not fixed
    K_fixed,            // Number of components fixed (i.e. at least one member has an observed label)
    L_inds,             // indices over views
    N_inds;             // indices over items
  
  arma::vec phis,
    mass, 
    complete_likelihood_vec, 
    N_ones, 
    N_log_factorial_vec;
  
  arma::umat 
    
    // Weight combination indices
    labels,
    
    // Various objects used to calculate MDI weights, normalising constant and phis
    comb_inds,
    phi_indicator,
    phi_ind_map,
    phi_indicator_t,
    
    // Class membership in each dataset
    N_k,
    
    // Indicator matrix for item n being an outlier in dataset l
    outliers,
    
    // Indicator matrix for item n being well-described by its component
    // in dataset l
    non_outliers,
    
    fixed;
  
  // The labels, weights in each dataset and the dataset correlation
  arma::mat w, phi_indicator_t_mat;
  
  // Cube of cluster members
  arma::ucube members;
  
  arma::field<arma::uvec> fixed_ind;
  
  // The data can have varying numbers of columns
  arma::field<arma::mat> X;
  
  // The collection of mixtures
  std::vector< std::unique_ptr<mixtureModel> > mixtures;
  // std::vector< mixtureModel > mixtures;
  
  mdiModelAlt(
    arma::field<arma::mat> _X,
    uvec _mixture_types,
    uvec _outlier_types,
    arma::uvec _K,
    arma::umat _labels,
    arma::umat _fixed
  ) ;
  
  // Destructor
  virtual ~mdiModelAlt() { };
  
  // === Components weights ====================================================
  
  // Calculate the rate parameter of the gamma distribution for the component 
  // weights
  double calcWeightRate(uword lstar, uword kstar);
  
  // Update the cluster weights
  void updateWeights();
  void updateWeightsViewL(uword l);
  
  void updateMassParameters();
  void updateMassParameterViewL(uword l);
  
  // === Phis ==================================================================
  
  // The rate for the phi coefficient between the lth and mth datasets.
  double calcPhiRate(uword lstar, uword mstar);
  
  // Calculate the log-weights of the mixture of Gammas the shape is sampled 
  // from
  arma::vec calculatePhiShapeMixtureWeights(int N_vw, double rate);
  
  // Sample the shape parameter of the phi posterior distribution from a mixture 
  // of Gammas
  double samplePhiShape(arma::uword l, arma::uword m, double rate);
  
  // Sample a new value for each phi parameter
  void updatePhis();
  
  // === Model parameters ======================================================
  
  // Updates the normalising constant for the posterior
  void updateNormalisingConst();
  
  // Sample a new value for the latent variable introduced to encourage nice 
  // posterior distributions
  void sampleStrategicLatentVariable();
  
  // Update the current allocations
  void updateAllocation();
  void updateAllocationViewL(uword l);
  
  // === Initialisation functions ==============================================
  
  // Initialise matrices relating to phis
  void initialisePhis();
  
  // This wrapper to declare the mixtures at the dataset level is kept as its
  // own separate function to make a semi-supervised class easier
  void initialiseMixtures();
  
  // Sample from the prior distribution of all parameters
  void sampleFromPriors();
  
  // Sample from the prior distribution of parameters at the view-specific level
  void sampleFromLocalPriors();
  
  // Sample from the prior distribution of parameters at the global level, i.e.
  // phis and gammas (arguably v and Z too)
  void sampleFromGlobalPriors();
  vec samplePhiPrior(uword n_phis);
  double sampleWeightPrior(uword l);
  vec sampleMassPrior();
  
  // Calculate the upweight due to matching labels across views (a function of 
  // phis)
  mat calculateUpweights(uword l);
  
  // Initialise the model, sampling from priors and calculating some initial 
  // objects
  void initialiseMDI();

  // === Label swapping functions ==============================================
  
  // This is used to consider possible label swaps
  double sampleLabel(arma::uword k, arma::vec K_inv_cum);
  
  // Calculate the score for a given labelling setup
  double calcScore(arma::uword l, arma::umat labels);
  
  // Check if labels should be swapped to improve correlation of clustering
  // across datasets via random sampling.
  arma::umat swapLabels(arma::uword l, arma::uword k, arma::uword k_prime);
  
  // Check if labels should be swapped to improve correlation of clustering
  // across datasets via random sampling.
  void updateLabels();
  void updateLabelsViewL(uword l);

};

#endif /* MDI_H */