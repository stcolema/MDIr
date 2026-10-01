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

//' @title Predict new items from saved MCMC draws
//' @description For each saved draw, loads the weights, phis and component
//' parameters and evaluates, for every new item, its marginal likelihood
//' (summed over all joint component assignments) and its posterior probability of
//' belonging to each component of each view. See `predictMDI()` in R.
//' @param X List of the data the model was fitted to (used to rebuild the
//' densities and the data-driven hyperparameters).
//' @param X_new List of the new data matrices, with the same columns as `X`.
//' @param K Number of components in each view.
//' @param mixture_types Integer density codes.
//' @param outlier_types Integer outlier component codes.
//' @param parameters For each view, a matrix with a row per draw holding the
//' flattened component parameters.
//' @param weights Cube (draws x K_max x L) of the component weights.
//' @param phis Matrix (draws x L(L - 1) / 2) of the phis.
//' @param outlier_weights Matrix (draws x L) of the outlier weights.
//' @param allocations Cube (draws x N x L) of the sampled labels of the fitted
//' items (0-based), used for the co-clustering probabilities.
//' @param coclustering Also return, for each view, the probability that each new
//' item shares a component with each fitted item.
//' @param prior Optional MDI-level prior vector (see `runMDI`).
//' @param density_prior Density-level prior options (see `runMDI`).
//' @return A list with `log_likelihood` (draws x new items), `class_probability`
//' (a K_max x new items matrix for each view, averaged over draws) and, if
//' requested, `coclustering` (a new items x N matrix for each view).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List predictNewItemsCpp(
    arma::field<arma::mat> X,
    arma::field<arma::mat> X_new,
    arma::uvec K,
    arma::uvec mixture_types,
    arma::uvec outlier_types,
    arma::field<arma::mat> parameters,
    arma::cube weights,
    arma::mat phis,
    arma::mat outlier_weights,
    arma::cube allocations,
    bool coclustering,
    arma::vec prior,
    arma::vec density_prior
) {
  std::unique_ptr<mdi> model = buildModel(X, K, mixture_types, outlier_types, prior, density_prior);
  const uword L = model->L, n_draws = weights.n_rows, n_new = X_new(0).n_rows;
  const uword N = allocations.n_cols;
  
  for(uword l = 0; l < L; l++) {
    if(X_new(l).n_rows != n_new) {
      Rcpp::stop("The new data must have the same number of items in every view.");
    }
    model->mixtures[l]->density_ptr->replaceData(X_new(l));
    model->mixtures[l]->outlierComponent_ptr->replaceData(X_new(l));
  }
  
  arma::mat log_lik(n_draws, n_new);
  arma::field<arma::mat> class_prob(L), cocluster(L);
  for(uword l = 0; l < L; l++) {
    class_prob(l).zeros(model->K_max, n_new);
    if(coclustering) {
      cocluster(l).zeros(n_new, N);
    }
  }
  
  arma::mat log_g;
  for(uword r = 0; r < n_draws; r++) {
    Rcpp::checkUserInterrupt();
    for(uword l = 0; l < L; l++) {
      for(uword k = 0; k < model->K_max; k++) {
        model->w(k, l) = weights(r, k, l);
      }
    }
    if(model->LC2 > 0 && L > 1) {
      model->phis = phis.row(r).t();
    }
    for(uword l = 0; l < L; l++) {
      auto& mixture = model->mixtures[l];
      mixture->density_ptr->setParameters(parameters(l).row(r).t());
      if(mixture->outlierComponent_ptr->active()) {
        mixture->outlierComponent_ptr->outlier_weight = outlier_weights(r, l);
        mixture->outlierComponent_ptr->non_outlier_weight = 1.0 - outlier_weights(r, l);
      }
    }
    model->refreshPartitionTables();
    const double log_Z = std::log(mdiPartitionSumFromC(model->w, model->K, model->partition_tables));
    
    for(uword j = 0; j < n_new; j++) {
      model->componentLogLikelihoods(j, false, log_g);
      double log_numerator = 0.0;
      const arma::mat probs = mdiClassProbabilities(model->w, model->K, model->partition_tables, log_g, log_numerator);
      log_lik(r, j) = log_numerator - log_Z;
      for(uword l = 0; l < L; l++) {
        class_prob(l).col(j) += probs.col(l) / (double) n_draws;
        if(coclustering) {
          for(uword i = 0; i < N; i++) {
            cocluster(l)(j, i) += probs((uword) std::llround(allocations(r, i, l)), l) / (double) n_draws;
          }
        }
      }
    }
  }
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("log_likelihood") = log_lik,
    Rcpp::Named("class_probability") = class_prob
  );
  if(coclustering) {
    out["coclustering"] = cocluster;
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

//' @title Test hook: log Metropolis-Hastings ratio of a label swap
//' @description Log acceptance ratio for exchanging two components of one view
//' (labels, weights and component parameters together; test hook).
//' @param labels N x L matrix of labels (0-based).
//' @param phi L x L matrix of phis.
//' @param w Weights (K_max x L).
//' @param K Components per view.
//' @param v Strategic latent variable.
//' @param lstar View (0-based).
//' @param k,kprime Components to exchange (0-based).
//' @return The log ratio.
//' @keywords internal
// [[Rcpp::export]]
double mdiSwapLogRatioCpp(arma::umat labels, arma::mat phi, arma::mat w, arma::uvec K, double v,
                          arma::uword lstar, arma::uword k, arma::uword kprime) {
  const std::vector<double> C = mdiConnectedSums(phi);
  double Z_swapped = 0.0;
  return mdiSwapLogRatio(labels, phi, w, K, C, v, lstar, k, kprime,
                         mdiPartitionSumFromC(w, K, C), Z_swapped);
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

//' @title Test hook: log marginal likelihood of an item
//' @description log Z(w * g) - log Z(w) for per-view component likelihoods
//' (test hook).
//' @param log_g K_max x L matrix of log component likelihoods (-Inf allowed).
//' @param w Weights (K_max x L). @param K Components per view. @param phi L x L matrix of phis.
//' @return The log marginal likelihood.
//' @keywords internal
// [[Rcpp::export]]
double mdiLogMarginalCpp(arma::mat log_g, arma::mat w, arma::uvec K, arma::mat phi) {
  const std::vector<double> C = mdiConnectedSums(phi);
  return mdiLogNumerator(w, K, C, log_g) - std::log(mdiPartitionSumFromC(w, K, C));
}

//' @title Test hook: class probabilities of an item
//' @description Posterior probability of each component of each view given an
//' item's data (test hook).
//' @param log_g K_max x L matrix of log component likelihoods (-Inf allowed).
//' @param w Weights (K_max x L). @param K Components per view. @param phi L x L matrix of phis.
//' @return A K_max x L matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::mat mdiClassProbabilitiesCpp(arma::mat log_g, arma::mat w, arma::uvec K, arma::mat phi) {
  double log_numerator = 0.0;
  return mdiClassProbabilities(w, K, mdiConnectedSums(phi), log_g, log_numerator);
}

//' @title Test hook: chain of draws from the collapsed conditional of a phi
//' @description Runs the slice-sampling update of phi(l, m) repeatedly with the
//' weights and the agreement count held fixed (test hook).
//' @param w Weights (K_max x L). @param K Components per view. @param phi L x L matrix of phis.
//' @param l,m Views (0-based).
//' @param N Number of items. @param N_lm Number of items with the same component in views l and m.
//' @param shape,rate Gamma prior on phi(l, m).
//' @param n_draws Number of updates.
//' @return The draws.
//' @keywords internal
// [[Rcpp::export]]
arma::vec phiSliceChainCpp(arma::mat w, arma::uvec K, arma::mat phi, arma::uword l, arma::uword m,
                           double N, double N_lm, double shape, double rate, arma::uword n_draws) {
  const double B = mdiPhiRate(w, K, phi, l, m);
  arma::mat phi_zero = phi;
  phi_zero(l, m) = 0.0;
  phi_zero(m, l) = 0.0;
  const double A = mdiPartitionSum(w, K, phi_zero);
  double current = phi(l, m);
  arma::vec out(n_draws);
  for(uword i = 0; i < n_draws; i++) {
    current = mdiSamplePhiSlice(current, N_lm, N, A, B, shape, rate);
    out(i) = current;
  }
  return out;
}

//' @title Test hook: log collapsed conditional of a phi
//' @description Unnormalised log density of phi(l, m) with the strategic latent
//' variable integrated out, as a function of phi (test hook).
//' @param phi Values of phi(l, m). @param w Weights. @param K Components per view.
//' @param phi_matrix L x L matrix of phis (the entry (l, m) is ignored).
//' @param l,m Views (0-based).
//' @param N Number of items. @param N_lm Agreement count.
//' @param shape,rate Gamma prior.
//' @return The log density (without the Jacobian of log phi).
//' @keywords internal
// [[Rcpp::export]]
arma::vec phiConditionalLogDensityCpp(arma::vec phi, arma::mat w, arma::uvec K, arma::mat phi_matrix,
                                      arma::uword l, arma::uword m, double N, double N_lm,
                                      double shape, double rate) {
  const double B = mdiPhiRate(w, K, phi_matrix, l, m);
  arma::mat phi_zero = phi_matrix;
  phi_zero(l, m) = 0.0;
  phi_zero(m, l) = 0.0;
  const double A = mdiPartitionSum(w, K, phi_zero);
  arma::vec out(phi.n_elem);
  for(uword i = 0; i < phi.n_elem; i++) {
    out(i) = mdiLogPhiConditional(std::log(phi(i)), N_lm, N, A, B, shape, rate) - std::log(phi(i));
  }
  return out;
}
