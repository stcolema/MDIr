// predictive.cpp
// =============================================================================
// Prior and posterior predictive simulation, and test hooks for the exact
// normalising-constant calculations.
# include "mdi.h"
# include "mdiPartition.h"
# include "mvt.h"

using namespace Rcpp ;
using namespace arma ;

namespace {

// A model object used only for its priors and simulators: labels are zero and 
// nothing is fixed.
std::unique_ptr<mdi> buildModel(
    const arma::field<arma::mat>& X,
    const arma::uvec& K,
    const arma::uvec& mixture_types,
    const arma::uvec& outlier_types,
    const arma::vec& prior,
    const arma::vec& density_prior
) {
  const uword L = X.n_elem, N = X(0).n_rows;
  arma::umat labels(N, L, arma::fill::zeros), fixed(N, L, arma::fill::zeros);
  return std::unique_ptr<mdi>(new mdi(X, mixture_types, outlier_types, K, labels, fixed, prior, density_prior));
}

}

//' @title Data-driven prior hyperparameters of a density
//' @description Constructs the density for the data `X` and returns the
//' hyperparameters of its priors (which for most densities are set from `X`).
//' @param X Data matrix. Non-finite entries are treated as missing.
//' @param K Number of components.
//' @param mixture_type Integer density code (0 = G, 1 = MVN, 2 = C, 3 = GP).
//' @param density_prior Density-level prior options (see `runMDI`).
//' @return A named list of hyperparameters.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List densityHyperparameters(arma::mat X, arma::uword K, arma::uword mixture_type, arma::vec density_prior) {
  arma::uvec labels(X.n_rows, arma::fill::zeros);
  densityFactory::densityType val = static_cast<densityFactory::densityType>(mixture_type);
  std::unique_ptr<density> d = densityFactory::createDensity(val, K, labels, X, density_prior);
  return d->hyperparameterList();
}

//' @title Simulate datasets from the prior predictive distribution
//' @description Draws the MDI parameters (masses, phis, weights), the labels 
//' and the component parameters from their priors, then simulates every view.
//' @param X List of data matrices; used only to set data-driven hyperparameters
//' and the dimensions. Non-finite entries are ignored.
//' @param K Number of components in each view.
//' @param mixture_types Integer density codes.
//' @param outlier_types Integer outlier component codes.
//' @param n_datasets Number of datasets to simulate.
//' @param prior Optional MDI-level prior vector (see `runMDI`).
//' @param density_prior Density-level prior options (see `runMDI`).
//' @return A list with one entry per dataset holding `data` (list of matrices), 
//' `labels`, `outliers`, `mass`, `phis` and `weights`.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List simulatePriorPredictiveCpp(
    arma::field<arma::mat> X,
    arma::uvec K,
    arma::uvec mixture_types,
    arma::uvec outlier_types,
    arma::uword n_datasets,
    arma::vec prior,
    arma::vec density_prior
) {
  std::unique_ptr<mdi> model = buildModel(X, K, mixture_types, outlier_types, prior, density_prior);
  const uword L = model->L, N = model->N;
  
  Rcpp::List out(n_datasets);
  for(uword d = 0; d < n_datasets; d++) {
    Rcpp::checkUserInterrupt();
    model->sampleFromGlobalPriors();
    model->sampleFromLocalPriors();
    for(uword l = 0; l < L; l++) {
      model->mixtures[l]->outlierComponent_ptr->sampleFromPrior();
    }
    
    const arma::umat labels = model->samplePriorLabels(N);
    arma::umat outliers(N, L, arma::fill::zeros);
    Rcpp::List data(L);
    
    for(uword l = 0; l < L; l++) {
      auto& mixture = model->mixtures[l];
      const bool has_outliers = mixture->outlierComponent_ptr->active();
      arma::mat sim(N, X(l).n_cols);
      for(uword n = 0; n < N; n++) {
        const bool outlier = has_outliers && (randu() < mixture->outlierComponent_ptr->outlier_weight);
        outliers(n, l) = outlier ? 1 : 0;
        sim.row(n) = (outlier ? mixture->outlierComponent_ptr->simulate() 
                              : mixture->density_ptr->simulate(labels(n, l))).t();
      }
      data[l] = sim;
    }
    out[d] = Rcpp::List::create(
      Rcpp::Named("data") = data,
      Rcpp::Named("labels") = labels,
      Rcpp::Named("outliers") = outliers,
      Rcpp::Named("mass") = model->mass,
      Rcpp::Named("phis") = model->phis,
      Rcpp::Named("weights") = model->w
    );
  }
  return out;
}

//' @title Simulate replicate datasets from the posterior predictive distribution
//' @description For each saved MCMC draw, loads that draw's component
//' parameters and simulates a replicate of every item from its sampled 
//' component (or the outlier distribution if it was sampled as an outlier).
//' @param X List of the observed data matrices (used to rebuild the densities).
//' @param K Number of components in each view.
//' @param mixture_types Integer density codes.
//' @param outlier_types Integer outlier component codes.
//' @param parameters For each view, a matrix with a row per draw holding the 
//' flattened component parameters.
//' @param allocations Cube (draws x N x L) of sampled labels (as doubles).
//' @param outliers Cube (draws x N x L) of sampled outlier indicators.
//' @param prior Optional MDI-level prior vector (see `runMDI`).
//' @param density_prior Density-level prior options (see `runMDI`).
//' @return A list with one entry per view: a cube (draws x N x P) of replicates.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List simulatePosteriorPredictiveCpp(
    arma::field<arma::mat> X,
    arma::uvec K,
    arma::uvec mixture_types,
    arma::uvec outlier_types,
    arma::field<arma::mat> parameters,
    arma::cube allocations,
    arma::cube outliers,
    arma::vec prior,
    arma::vec density_prior
) {
  std::unique_ptr<mdi> model = buildModel(X, K, mixture_types, outlier_types, prior, density_prior);
  const uword L = model->L, N = model->N, n_draws = allocations.n_rows;
  
  Rcpp::List out(L);
  for(uword l = 0; l < L; l++) {
    auto& mixture = model->mixtures[l];
    const bool has_outliers = mixture->outlierComponent_ptr->active();
    arma::cube replicates(n_draws, N, X(l).n_cols);
    for(uword r = 0; r < n_draws; r++) {
      Rcpp::checkUserInterrupt();
      mixture->density_ptr->setParameters(parameters(l).row(r).t());
      for(uword n = 0; n < N; n++) {
        const bool outlier = has_outliers && outliers(r, n, l) > 0.5;
        const arma::vec x = outlier ? mixture->outlierComponent_ptr->simulate()
                                    : mixture->density_ptr->simulate((uword) std::llround(allocations(r, n, l)));
        for(uword p = 0; p < x.n_elem; p++) {
          replicates(r, n, p) = x(p);
        }
      }
    }
    out[l] = replicates;
  }
  return out;
}

// Test hooks for the exact normalising-constant calculations ------------------

//' @title Test hook: MDI normalising constant
//' @description Exact normalising constant Z of the MDI model (test hook).
//' @param w Weights (K_max x L). @param K Components per view. @param phi L x L matrix of phis.
//' @return Z.
//' @keywords internal
// [[Rcpp::export]]
double mdiNormalisingConstantCpp(arma::mat w, arma::uvec K, arma::mat phi) {
  return mdiPartitionSum(w, K, phi);
}

//' @title Test hook: rate of a weight conditional
//' @description dZ/dw for one weight (test hook).
//' @param w Weights. @param K Components per view. @param phi L x L matrix of phis.
//' @param lstar View (0-based). @param kstar Component (0-based).
//' @return The derivative.
//' @keywords internal
// [[Rcpp::export]]
double mdiWeightRateCpp(arma::mat w, arma::uvec K, arma::mat phi, arma::uword lstar, arma::uword kstar) {
  return mdiWeightRate(w, K, phi, lstar, kstar);
}

//' @title Test hook: rates of all the weight conditionals of a view
//' @description dZ/dw for every weight of one view from a single pass (test hook).
//' @param w Weights (K_max x L).
//' @param K Components per view.
//' @param phi L x L matrix of phis.
//' @param lstar View (0-based).
//' @return The derivatives for components 0, ..., K(lstar) - 1.
//' @keywords internal
// [[Rcpp::export]]
std::vector<double> mdiWeightRatesCpp(arma::mat w, arma::uvec K, arma::mat phi, arma::uword lstar) {
  return arma::conv_to< std::vector<double> >::from(mdiWeightRates(w, K, mdiConnectedSums(phi), lstar));
}

//' @title Test hook: rate of a phi conditional
//' @description dZ/dphi for one pair of views (test hook).
//' @param w Weights. @param K Components per view. @param phi L x L matrix of phis.
//' @param l,m Views (0-based).
//' @return The derivative.
//' @keywords internal
// [[Rcpp::export]]
double mdiPhiRateCpp(arma::mat w, arma::uvec K, arma::mat phi, arma::uword l, arma::uword m) {
  return mdiPhiRate(w, K, phi, l, m);
}

//' @title Test hook: impute from the outlier distribution
//' @description Draws `n_rep` observations from the multivariate t outlier 
//' distribution defined by `X`, hides the entries `missing_cols` (1-based) and 
//' re-imputes them from their exact conditional. If the imputation is exact, 
//' the returned matrix has the same distribution as the direct draws in `direct`.
//' @param X Data defining the t location and scale (finite entries).
//' @param missing_cols 1-based columns to hide and impute.
//' @param n_rep Number of replicates.
//' @return List with `direct` and `imputed` matrices (n_rep x P).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mvtImputationCheckCpp(arma::mat X, arma::uvec missing_cols, arma::uword n_rep) {
  const uword P = X.n_cols;
  arma::field<arma::uvec> miss(X.n_rows), obs(X.n_rows);
  arma::uvec miss_idx = missing_cols - 1;
  arma::uvec obs_idx;
  {
    arma::uvec all = arma::regspace<arma::uvec>(0, P - 1);
    std::vector<uword> keep;
    for(uword p = 0; p < P; p++) {
      if(arma::accu(miss_idx == p) == 0) keep.push_back(p);
    }
    obs_idx = arma::uvec(keep);
  }
  for(uword n = 0; n < X.n_rows; n++) {
    miss(n) = miss_idx;
    obs(n) = obs_idx;
  }
  
  arma::uvec fixed(X.n_rows, arma::fill::zeros);
  mvt component(fixed, X, &miss, &obs);
  
  // A one-row copy of the component whose "data" is replaced draw by draw
  arma::mat direct(n_rep, P), imputed(n_rep, P);
  for(uword r = 0; r < n_rep; r++) {
    const arma::vec x = component.simulate();
    direct.row(r) = x.t();
    component.X = x.t();
    const arma::vec filled = component.sampleMissingValues(0);
    arma::vec x_new = x;
    for(uword i = 0; i < miss_idx.n_elem; i++) {
      x_new(miss_idx(i)) = filled(i);
    }
    imputed.row(r) = x_new.t();
  }
  return Rcpp::List::create(Rcpp::Named("direct") = direct, Rcpp::Named("imputed") = imputed);
}

//' @title Test hook: inverse-gamma calibration
//' @description Inverse-gamma parameters with `tail` probability below `lower` and
//' above `upper` (used for the Gaussian process length scale prior).
//' @param lower,upper Bounds.
//' @param tail Tail probability.
//' @return c(shape, rate).
//' @keywords internal
// [[Rcpp::export]]
arma::vec calibrateInverseGammaCpp(double lower, double upper, double tail) {
  return calibrateInverseGamma(lower, upper, tail);
}

//' @title Test hook: population update of the GP hyperparameters
//' @description Runs the update of the mean and sd of a population of log
//' hyperparameters `y` repeatedly with `y` held fixed, so the draws can be 
//' compared with the analytic posterior of (m, s) given `y`.
//' @param y Log hyperparameters of the occupied components.
//' @param centre Prior mean of the population mean.
//' @param center_sd Prior sd of the population mean.
//' @param pool_sd_scale Scale of the half-normal prior on the population sd.
//' @param n_iter Number of updates.
//' @return A matrix with columns m and s.
//' @keywords internal
// [[Rcpp::export]]
arma::mat gpPopulationCheckCpp(arma::vec y, double centre, double center_sd, 
                               double pool_sd_scale, arma::uword n_iter) {
  arma::mat X(6, 3, arma::fill::randn);
  arma::uvec labels(6, arma::fill::zeros);
  gp g(1, labels, X, arma::vec({2.0, 1.0, 1.0, center_sd, pool_sd_scale}));
  g.log_variance_centre = centre;
  double m = centre, s = 1.0;
  arma::mat out(n_iter, 2);
  for(uword it = 0; it < n_iter; it++) {
    g.updatePopulation(y, m, s);
    out(it, 0) = m;
    out(it, 1) = s;
  }
  return out;
}
