// mdi.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "mdi.h"

using namespace arma ;

// =============================================================================
// MDI class

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;


mdi::mdi(
  arma::field<arma::mat> _X,
  uvec _mixture_types,
  uvec _outlier_types,
  arma::uvec _K,
  arma::umat _labels,
  arma::umat _fixed,
  arma::vec _prior,
  arma::vec _density_prior
) {

  setPrior(_prior);
  density_prior = _density_prior;
  mixture_types = _mixture_types;
  outlier_types = _outlier_types;

  // The number of datasets modelled
  L = _X.n_elem;
  if(L < 1 || L > 16) {
    Rcpp::stop("The number of views must be between 1 and 16.");
  }

  // The number of pairwise combinations
  if(L > 1) {
    LC2 = L * (L - 1) / 2;
  }

  // All views must describe the same items
  N = _X(0).n_rows;
  for(uword l = 0; l < L; l++) {
    if(_X(l).n_rows != N) {
      Rcpp::stop("Datasets not matching in number of rows.");
    }
  }

  // The number of components modelled in each dataset
  K = _K;
  if(K.n_elem != L || any(K < 1)) {
    Rcpp::stop("K must hold a positive number of components for each view.");
  }
  K_max = max(K);

  N_k.set_size(K_max, L);
  N_k.zeros();

  mass.set_size(L);
  mass.zeros();
  mass_acceptance_count.set_size(L);
  mass_acceptance_count.zeros();

  // A weight vector for each dataset. For ease of manipulation K_max is used
  // for every view to avoid ragged fields of vectors.
  w.set_size(K_max, L);
  w.zeros();

  phis.set_size(LC2);
  phis.zeros();

  // Map between a dataset pair and the index of its phi
  phi_map.set_size(L, L);
  phi_map.zeros();
  uword col_ind = 0;
  for(uword l = 0; l + 1 < L; l++) {
    for(uword m = l + 1; m < L; m++) {
      phi_map(m, l) = col_ind;
      phi_map(l, m) = col_ind;
      col_ind++;
    }
  }

  labels = _labels;
  X = _X;

  members.set_size(N, K_max, L);
  members.zeros();

  non_outliers.set_size(N, L);
  non_outliers.ones();
  outliers.set_size(N, L);
  outliers.zeros();

  // The observed labels
  fixed = _fixed;

  // The components holding an observed label are fixed; the others are free. The
  // observed labels never change, so this is computed once. Nothing is assumed about
  // which indices the observed classes occupy.
  free_components.assign(L, arma::uvec());
  for(uword l = 0; l < L; l++){
    arma::uvec has_observed(K(l), arma::fill::zeros);
    for(uword n = 0; n < N; n++) {
      if(fixed(n, l) == 1) {
        if(labels(n, l) >= K(l)) {
          Rcpp::stop("An observed label lies outside 0, ..., K - 1 in view %d.", (int) l + 1);
        }
        has_observed(labels(n, l)) = 1;
      }
    }
    free_components[l] = arma::find(has_observed == 0);
  }

  complete_likelihood_vec = zeros< vec >(L);
  observed_likelihood_vec = zeros< vec >(L);

  initialiseMDI();
};

void mdi::setPrior(const arma::vec& prior) {
  if(prior.n_elem == 0) {
    return;
  }
  if(prior.n_elem != 5 || !all(prior > 0.0) || !prior.is_finite()) {
    Rcpp::stop("The prior must hold five positive values: mass shape and rate, weight rate, phi shape and rate.");
  }
  mass_shape_prior = prior(0);
  mass_rate_prior = prior(1);
  w_rate_prior = prior(2);
  phi_shape_prior = prior(3);
  phi_rate_prior = prior(4);
}

arma::mat mdi::phiMatrix() const {
  arma::mat phi_mat(L, L, arma::fill::zeros);
  for(uword l = 0; l + 1 < L; l++) {
    for(uword m = l + 1; m < L; m++) {
      phi_mat(l, m) = phis(phi_map(l, m));
      phi_mat(m, l) = phi_mat(l, m);
    }
  }
  return phi_mat;
}

void mdi::initialiseMixtures() {
  mixtures.reserve(L);
  for(uword l = 0; l < L; l++) {
    mixtures.push_back(
      std::unique_ptr<mixtureModel>(new mixtureModel(
        mixture_types(l),
        outlier_types(l),
        K(l),
        labels.col(l),
        fixed.col(l),
        X(l),
        density_prior
      ))
    );
    non_outliers.col(l) = mixtures[l]->non_outliers;
    outliers.col(l) = mixtures[l]->outliers;
  }
};

// === Normalising constant, weights and phis ==================================

double mdi::calcPhiRate(uword l, uword m) const {
  return v * mdiPhiRate(w, K, phiMatrix(), l, m);
}

void mdi::refreshPartitionTables() {
  if(partition_tables.empty()
       || partition_tables_phis.n_elem != phis.n_elem
       || !all(partition_tables_phis == phis)) {
    partition_tables = mdiConnectedSums(phiMatrix());
    partition_tables_phis = phis;
  }
}

void mdi::updateNormalisingConstant() {
  refreshPartitionTables();
  Z = mdiPartitionSumFromC(w, K, partition_tables);
}

void mdi::sampleStrategicLatentVariable() {
  v = rGamma((double) N, Z);
};

void mdi::updateWeightsViewL(uword l) {

  // The phi matrix and the other views' weights are fixed while the weights of
  // view l are updated, and w(k, l) enters Z linearly, so each weight has a
  // Gamma full conditional given the others. The rate of w(k, l) does not
  // involve view l, so all K(l) rates come from a single pass.
  refreshPartitionTables();
  const arma::vec rates = mdiWeightRates(w, K, partition_tables, l);

  for(uword k = 0; k < K(l); k++) {
    double posterior_shape = (mass(l) / (double) K(l)) + (double) N_k(k, l);
    double posterior_rate = w_rate_prior + v * rates(k);
    w(k, l) = rGamma(posterior_shape, posterior_rate);
  }
}

void mdi::updateWeights() {
  for(uword l = 0; l < L; l++) {
    updateWeightsViewL(l);
  }
};

arma::vec mdi::calculatePhiShapeMixtureWeights(uword N_lm, double rate) const {

  // Marginalising the binomial expansion of prod_n (1 + phi 1[c_nl = c_nm])
  // gives a mixture over r = 0..N_lm of Gamma(phi_shape_prior + r,
  // phi_rate_prior + rate) with these weights.
  vec log_weights(N_lm + 1);
  const double log_rate = std::log(rate + phi_rate_prior);
  for(uword r = 0; r <= N_lm; r++) {
    log_weights(r) = logChoose((double) N_lm, (double) r)
      + std::lgamma((double) r + phi_shape_prior)
      - ((double) r + phi_shape_prior) * log_rate;
  }
  return log_weights;
};

uword mdi::samplePhiShape(uword N_lm, double rate) const {
  vec log_weights = calculatePhiShapeMixtureWeights(N_lm, rate);
  vec weights = exp(log_weights - max(log_weights));
  weights /= accu(weights);
  return sampleCategorical(weights);
}

void mdi::updatePhis() {
  if(L == 1) {
    return;
  }
  for(uword l = 0; l + 1 < L; l++) {
    for(uword m = l + 1; m < L; m++) {
      double rate = calcPhiRate(l, m);
      uword N_lm = accu(labels.col(l) == labels.col(m));
      uword r = samplePhiShape(N_lm, rate);
      phis(phi_map(l, m)) = rGamma(phi_shape_prior + (double) r, phi_rate_prior + rate);
    }
  }
};

void mdi::updatePhisSlice() {
  if(L == 1) {
    return;
  }
  for(uword l = 0; l + 1 < L; l++) {
    for(uword m = l + 1; m < L; m++) {
      const uword idx = phi_map(l, m);
      const double N_lm = (double) accu(labels.col(l) == labels.col(m));

      // Z is linear in phi(l, m): Z = A + B phi
      arma::mat phi_mat = phiMatrix();
      const double B = mdiPhiRate(w, K, phi_mat, l, m);
      phi_mat(l, m) = 0.0;
      phi_mat(m, l) = 0.0;
      const double A = mdiPartitionSum(w, K, phi_mat);

      phis(idx) = mdiSamplePhiSlice(phis(idx), N_lm, (double) N, A, B, phi_shape_prior, phi_rate_prior);
    }
  }
}

// === Likelihood ==============================================================

void mdi::componentLogLikelihoods(uword n, bool use_fixed, arma::mat& log_g) {
  log_g.set_size(K_max, L);
  log_g.fill(-arma::datum::inf);
  for(uword l = 0; l < L; l++) {
    auto& mixture = mixtures[l];
    const arma::vec ll = mixture->density_ptr->itemLogLikelihood(n);
    if(use_fixed && fixed(n, l) == 1) {
      // An observed label: the item belongs to that component, never to the outlier
      log_g(labels(n, l), l) = ll(labels(n, l));
    } else if(mixture->outlierComponent_ptr->active()) {
      const double log_w_non = std::log(mixture->outlierComponent_ptr->non_outlier_weight);
      const double log_w_out = std::log(mixture->outlierComponent_ptr->outlier_weight);
      const double ll_out = log_w_out + mixture->outlierComponent_ptr->outlier_likelihood(n);
      for(uword k = 0; k < K(l); k++) {
        log_g(k, l) = logSumExp(arma::vec({log_w_non + ll(k), ll_out}));
      }
    } else {
      for(uword k = 0; k < K(l); k++) {
        log_g(k, l) = ll(k);
      }
    }
  }
}

arma::vec mdi::pointwiseLogLikelihood() {
  refreshPartitionTables();
  const double log_Z = std::log(mdiPartitionSumFromC(w, K, partition_tables));

  arma::vec out(N);
  arma::mat log_g;
  for(uword n = 0; n < N; n++) {
    componentLogLikelihoods(n, true, log_g);
    out(n) = mdiLogNumerator(w, K, partition_tables, log_g) - log_Z;
  }
  return out;
}

// === Tempering ===============================================================

void mdi::setBeta(double beta_new) {
  for(uword l = 0; l < L; l++) {
    mixtures[l]->setBeta(beta_new);
  }
  beta = beta_new;
}

double mdi::dataLogLikelihood() {
  double total = 0.0;
  for(uword l = 0; l < L; l++) {
    const auto& density = mixtures[l]->density_ptr;
    for(uword n = 0; n < N; n++) {
      total += density->logLikelihood(n, labels(n, l));
    }
  }
  return total;
}

mdiState mdi::saveState() const {
  mdiState s;
  s.labels = labels;
  s.w = w;
  s.phis = phis;
  s.mass = mass;
  s.v = v;
  s.Z = Z;
  for(uword l = 0; l < L; l++) {
    s.theta.push_back(mixtures[l]->density_ptr->parameters());
    s.pooled.push_back(mixtures[l]->density_ptr->pooledHyperparameters());
  }
  return s;
}

void mdi::loadState(const mdiState& s) {
  labels = s.labels;
  w = s.w;
  phis = s.phis;
  mass = s.mass;
  v = s.v;
  Z = s.Z;
  partition_tables.clear();
  partition_tables_phis.reset();
  for(uword l = 0; l < L; l++) {
    auto& density = mixtures[l]->density_ptr;
    density->setParameters(s.theta[l]);
    if(s.pooled[l].n_elem > 0) {
      density->setPooledHyperparameters(s.pooled[l]);
    }
    mixtures[l]->labels = labels.col(l);
    density->labels = labels.col(l);
    refreshMembersViewL(l);
  }
}

void mdi::initialiseFromPrior() {
  sampleFromGlobalPriors();
  sampleFromLocalPriors();
  labels = samplePriorLabels(N);
  for(uword l = 0; l < L; l++) {
    mixtures[l]->labels = labels.col(l);
    mixtures[l]->density_ptr->labels = labels.col(l);
    refreshMembersViewL(l);
  }
  Z_start = Z;
}

// === Priors ==================================================================

void mdi::sampleFromPriors() {
  sampleFromGlobalPriors();
  sampleFromLocalPriors();
};

void mdi::sampleFromLocalPriors() {
  for(uword l = 0; l < L; l++) {
    mixtures[l]->sampleFromPriors();
  }
};

vec mdi::samplePhiPrior(uword n_phis) {
  return rGamma(n_phis, phi_shape_prior , phi_rate_prior);
};

double mdi::sampleWeightPrior(uword l) {
  return rGamma(mass(l) / (double) K(l) , w_rate_prior);
}

vec mdi::sampleMassPrior() {
  return rGamma(L, mass_shape_prior, mass_rate_prior);
}

void mdi::sampleFromGlobalPriors() {
  mass = sampleMassPrior();
  if(L > 1) {
    phis = samplePhiPrior(LC2);
  } else {
    phis.zeros();
  }
  w.zeros();
  for(uword l = 0; l < L; l++) {
    for(uword k = 0; k < K(l); k++) {
      w(k, l) = sampleWeightPrior(l);
    }
  }
  updateNormalisingConstant();
  sampleStrategicLatentVariable();
};

// === Mass ====================================================================

void mdi::updateMassParameters() {
  for(uword l = 0; l < L; l++) {
    updateMassParameterViewL(l);
  }
}

void mdi::updateMassParameterViewL(uword lstar) {

  const double current_mass = mass(lstar);
  const vec current_weights = w(span(0, K(lstar) - 1), lstar);

  // Log-scale random walk. The target is a density in the mass itself, so the
  // Jacobian of the log transformation enters the acceptance ratio.
  const double proposed_mass = current_mass * std::exp(mass_proposal_sd * randn());
  if(!(proposed_mass > 0.0) || !std::isfinite(proposed_mass)) {
    return;
  }

  auto log_target = [&](double m) {
    return gammaLogLikelihood(current_weights, m / (double) K(lstar), w_rate_prior)
      + gammaLogLikelihood(m, mass_shape_prior, mass_rate_prior);
  };

  const double log_ratio = log_target(proposed_mass) - log_target(current_mass)
    + std::log(proposed_mass) - std::log(current_mass);

  if(std::log(randu()) < log_ratio) {
    mass(lstar) = proposed_mass;
    mass_acceptance_count(lstar)++;
  }
};

// === Allocations =============================================================

mat mdi::calculateUpweights(uword lstar) const {
  mat log_upweights(K(lstar), N, arma::fill::zeros);
  for(uword m = 0; m < L; m++) {
    if(m != lstar) {
      const double log_up = std::log1p(phis(phi_map(m, lstar)));
      for(uword n = 0; n < N; n++) {
        const uword k = labels(n, m);
        if(k < K(lstar)) {
          log_upweights(k, n) += log_up;
        }
      }
    }
  }
  return log_upweights;
};

void mdi::initialiseDatasetL(uword l) {
  vec log_weights = log(w(span(0, K(l) - 1), l));
  mat log_upweights = calculateUpweights(l);
  mixtures[l]->sampleAllMissingValues();
  mixtures[l]->initialiseMixture(log_weights, log_upweights);
  labels.col(l) = mixtures[l]->labels;
  non_outliers.col(l) = mixtures[l]->non_outliers;
  outliers.col(l) = mixtures[l]->outliers;
  // Likelihoods of the initial state
  complete_likelihood_vec(l) = mixtures[l]->complete_likelihood;
  observed_likelihood_vec(l) = mixtures[l]->observed_likelihood;
  complete_likelihood = accu(complete_likelihood_vec);
  observed_likelihood = accu(observed_likelihood_vec);
  refreshMembersViewL(l);
}

void mdi::sweep(uword iteration) {
  updateNormalisingConstant();
  Z_start = Z;
  if(phi_slice) {
    // The collapsed update moves phi, so Z must be refreshed before v is drawn
    updatePhisSlice();
    updateNormalisingConstant();
  }
  sampleStrategicLatentVariable();
  updateMassParameters();
  updateWeights();
  if(!phi_slice) {
    updatePhis();
  }
  for(uword l = 0; l < L; l++) {
    mixtures[l]->sampleParameters();
  }
  updateAllocation();
  if(split_merge_moves > 0) {
    updateSplitMerge();
  }
  if((iteration + 1) % 10 == 0) {
    updateLabels();
  }
}

arma::umat mdi::samplePriorLabels(uword n_items) const {
  double n_combinations = 1.0;
  for(uword l = 0; l < L; l++) {
    n_combinations *= (double) K(l);
  }
  if(n_combinations > 5e6) {
    Rcpp::stop("Too many joint component combinations (%g) to enumerate.", n_combinations);
  }
  const uword n_comb = (uword) n_combinations;
  
  arma::umat combinations(n_comb, L);
  arma::vec log_prob(n_comb);
  const arma::mat phi_mat = phiMatrix();
  for(uword i = 0; i < n_comb; i++) {
    uword remainder = i;
    for(uword l = 0; l < L; l++) {
      combinations(i, l) = remainder % K(l);
      remainder /= K(l);
    }
    double lp = 0.0;
    for(uword l = 0; l < L; l++) {
      lp += std::log(w(combinations(i, l), l));
    }
    for(uword l = 0; l + 1 < L; l++) {
      for(uword m = l + 1; m < L; m++) {
        if(combinations(i, l) == combinations(i, m)) {
          lp += std::log1p(phi_mat(l, m));
        }
      }
    }
    log_prob(i) = lp;
  }
  const arma::vec prob = arma::exp(log_prob - logSumExp(log_prob));
  
  arma::umat out(n_items, L);
  for(uword n = 0; n < n_items; n++) {
    out.row(n) = combinations.row(sampleCategorical(prob));
  }
  return out;
}

void mdi::initialiseMDI() {
  initialiseMixtures();
  sampleFromPriors();
  Z_start = Z;
  for(uword l = 0; l < L; l++) {
    initialiseDatasetL(l);
  }
};

void mdi::refreshMembersViewL(uword l) {
  members.slice(l).zeros();
  N_k.col(l).zeros();
  for(uword n = 0; n < N; n++) {
    members(n, labels(n, l), l) = 1;
    N_k(labels(n, l), l)++;
  }
  mixtures[l]->members = members.slice(l).cols(0, K(l) - 1);
  mixtures[l]->N_k = N_k(span(0, K(l) - 1), l);
}

void mdi::updateAllocationViewL(uword l) {
  vec log_weights = log(w(span(0, K(l) - 1), l));
  mat log_upweights = calculateUpweights(l);

  // Updates the labels, outlier indicators and (after the labels) any missing
  // values within the mixture using the MDI level weights and phis
  mixtures[l]->updateAllocation(log_weights, log_upweights);

  labels.col(l) = mixtures[l]->labels;
  non_outliers.col(l) = mixtures[l]->non_outliers;
  outliers.col(l) = mixtures[l]->outliers;
  complete_likelihood_vec(l) = mixtures[l]->complete_likelihood;
  observed_likelihood_vec(l) = mixtures[l]->observed_likelihood;
  refreshMembersViewL(l);
}

void mdi::updateAllocation() {
  for(uword l = 0; l < L; l++) {
    updateAllocationViewL(l);
  }
  complete_likelihood = accu(complete_likelihood_vec);
  observed_likelihood = accu(observed_likelihood_vec);
};

// === Label swapping ==========================================================

void mdi::updateLabels() {
  if(L == 1) {
    return;
  }
  for(uword l = 0; l < L; l++) {
    updateLabelsViewL(l);
  }
};

void mdi::updateLabelsViewL(uword lstar) {

  const arma::uvec& free_k = free_components[lstar];
  const uword n_free = free_k.n_elem;
  if(n_free < 2) {
    return;
  }

  refreshPartitionTables();
  const arma::mat phi_mat = phiMatrix();
  double Z_current = mdiPartitionSumFromC(w, K, partition_tables);

  for(uword i = 0; i < n_free; i++) {

    // Choose another free component uniformly. The proposal is symmetric.
    const uword k = free_k(i);
    uword j = std::min<uword>((uword) std::floor(randu() * (double) (n_free - 1)), n_free - 2);
    if(j >= i) {
      j++;
    }
    const uword k_prime = free_k(j);

    if(N_k(k, lstar) == 0 && N_k(k_prime, lstar) == 0) {
      continue;
    }

    // Exchanging (labels, weights, component parameters) of this view leaves the
    // prior over weights and parameters, the weight products and the data
    // likelihood unchanged. What changes is the phi coupling term and Z, which
    // depends on the weights through the alignment with the other views.
    double Z_swapped = 0.0;
    const double log_acceptance = mdiSwapLogRatio(
      labels, phi_mat, w, K, partition_tables, v, lstar, k, k_prime, Z_current, Z_swapped
    );

    if(std::log(randu()) < log_acceptance) {
      acceptance_count++;
      uvec loc_labs = labels.col(lstar);
      uvec in_k = find(loc_labs == k), in_k_prime = find(loc_labs == k_prime);
      loc_labs.elem(in_k).fill(k_prime);
      loc_labs.elem(in_k_prime).fill(k);
      labels.col(lstar) = loc_labs;
      std::swap(w(k, lstar), w(k_prime, lstar));
      Z_current = Z_swapped;
      mixtures[lstar]->labels = labels.col(lstar);
      mixtures[lstar]->swapComponents(k, k_prime);
      refreshMembersViewL(lstar);
    }
  }
};


// === Split-merge ==============================================================

void mdi::setSplitMerge(uword moves) {
  if(moves > 0) {
    for(uword l = 0; l < L; l++) {
      const auto& mix = mixtures[l];
      if(!mix->density_ptr->hasCollapsedMarginal()) {
        Rcpp::stop("The split-merge move needs 'G', 'MVN' or 'C' views (a collapsed marginal likelihood).");
      }
    }
  }
  split_merge_moves = moves;
}

void mdi::updateSplitMergeViewL(uword l) {
  if(K(l) < 2) {
    return;
  }
  const vec log_weights = log(w(span(0, K(l) - 1), l));
  const mat log_upweights = calculateUpweights(l);
  auto& mixture = mixtures[l];
  for(uword move = 0; move < split_merge_moves; move++) {
    // two distinct components, uniformly
    const uword a = std::min<uword>((uword) std::floor(randu() * (double) K(l)), K(l) - 1);
    uword b = std::min<uword>((uword) std::floor(randu() * (double) (K(l) - 1)), K(l) - 2);
    if(b >= a) {
      b++;
    }
    const bool accepted = mixture->splitMergeMove(a, b, log_weights, log_upweights);
    split_merge_attempts++;
    if(accepted) {
      split_merge_accepts++;
    }
    labels.col(l) = mixture->labels;
    non_outliers.col(l) = mixture->non_outliers;
    outliers.col(l) = mixture->outliers;
    refreshMembersViewL(l);
    // Redraw the parameters of the two components given the new labels (also after
    // a rejection: the move is a kernel on the labels with the parameters collapsed,
    // followed by a draw of the parameters from their conditional)
    mixture->density_ptr->labels = mixture->labels;
    mixture->density_ptr->sampleKthComponentParameters(a, mixture->members, mixture->non_outliers);
    mixture->density_ptr->sampleKthComponentParameters(b, mixture->members, mixture->non_outliers);
  }
  mixture->refreshLikelihoods(log_weights);
  complete_likelihood_vec(l) = mixture->complete_likelihood;
  observed_likelihood_vec(l) = mixture->observed_likelihood;
}

void mdi::updateSplitMerge() {
  for(uword l = 0; l < L; l++) {
    updateSplitMergeViewL(l);
  }
  complete_likelihood = accu(complete_likelihood_vec);
  observed_likelihood = accu(observed_likelihood_vec);
}
