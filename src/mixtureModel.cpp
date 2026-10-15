// mixtureModel.cpp
// =============================================================================
# include <RcppArmadillo.h>
# include "mixtureModel.h"

using namespace arma ;

// =============================================================================
// mixtureModel class

mixtureModel::mixtureModel(
  arma::uword _mixture_type,
  arma::uword _outlier_type,
  arma::uword _K,
  arma::uvec _labels,
  arma::uvec _fixed,
  arma::mat _X,
  arma::vec _density_prior) {
  
  mixture_type = _mixture_type;
  outlier_type = _outlier_type;
  K = _K;
  labels = _labels;
  
  N = _X.n_rows;
  P = _X.n_cols;
  
  if(labels.n_elem != N || _fixed.n_elem != N) {
    Rcpp::stop("Labels and fixed indicators must have one entry per item.");
  }
  if(labels.max() >= K) {
    Rcpp::stop("Initial labels must lie in 0, ..., K - 1.");
  }
  
  N_k = zeros<uvec>(K);
  complete_likelihood_vec = zeros< vec >(N);
  observed_likelihood_vec = zeros< vec >(N);
  members = zeros<umat>(N, K);
  alloc = zeros<mat>(N, K);
  
  fixed = _fixed;
  N_fixed = accu(fixed);
  
  // Observed labels have allocation probability one
  for(uword n = 0; n < N; n++) {
    if(fixed(n) == 1) {
      alloc(n, labels(n)) = 1.0;
    }
  }
  
  outliers = zeros<uvec>(N);
  non_outliers = ones<uvec>(N);
  
  // The density holds the data and the missing value patterns; the outlier
  // component refers to those patterns
  initialiseDensity(mixture_type, _X, _density_prior);
  n_param = density_ptr->n_param;
  initialiseOutlierComponent(outlier_type, _X);
};

void mixtureModel::initialiseDensity(arma::uword type, const arma::mat& X, const arma::vec& density_prior) {
  densityFactory::densityType val = static_cast<densityFactory::densityType>(type);
  density_ptr = densityFactory::createDensity(val, K, labels, X, density_prior);
};

void mixtureModel::initialiseOutlierComponent(arma::uword type, const arma::mat& X) {
  outlierComponentFactory::outlierType val = static_cast<outlierComponentFactory::outlierType>(type);
  outlierComponent_ptr = outlierComponentFactory::createOutlierComponent(
    val, fixed, X, &density_ptr->missing_indices, &density_ptr->observed_indices
  );
  outliers = outlierComponent_ptr->outliers;
  non_outliers = outlierComponent_ptr->non_outliers;
};

void mixtureModel::sampleFromPriors() {
  density_ptr->sampleFromPriors();
};

void mixtureModel::sampleParameters() {
  density_ptr->sampleParameters(members, non_outliers);
};

void mixtureModel::updateOutlierWeights() {
  non_outliers = 1 - outliers;

  // Only items without an observed label can be outliers, and only they carry a
  // factor of the outlier weight in the likelihood, so only they enter the Beta
  // update. (Observed items are never outliers; counting them as non-outliers would
  // bias the weight towards zero.)
  const arma::uvec is_free = 1 - fixed;
  outlierComponent_ptr->updateWeights(non_outliers % is_free, outliers % is_free);
}

void mixtureModel::updateAllocation(const arma::vec& log_weights, const arma::mat& log_upweights) {
  complete_likelihood = 0.0;
  observed_likelihood = 0.0;
  
  updateOutlierWeights();
  
  // The likelihood of each item uses its observed entries only, so the new 
  // (label, outlier) draw does not depend on the current imputations. The 
  // imputations are refreshed straight afterwards, conditional on the new 
  // allocation, before any parameter is updated.
  for(uword n = 0; n < N; n++) {
    updateItemAllocation(n, log_weights, log_upweights.col(n));
  }
  non_outliers = 1 - outliers;
  sampleAllMissingValues();
  
  observed_likelihood = accu(observed_likelihood_vec);
  complete_likelihood = accu(complete_likelihood_vec);
  
  { uvec occupied = unique(labels); K_occ = occupied.n_elem; }
}  

void mixtureModel::updateItemAllocation(
    uword n, 
    const arma::vec& log_weights, 
    const arma::vec& log_upweights
) {
  const bool has_outliers = outlierComponent_ptr->active();
  const vec ll = density_ptr->itemLogLikelihood(n);
  const vec log_component_weight = log_weights + log_upweights;
  
  // The allocation targets the tempered conditional, in which the component 
  // log-likelihood enters multiplied by beta. The recorded likelihoods below
  // stay untempered.
  const vec ll_sampling = (beta == 1.0) ? ll : vec(beta * ll);
  
  const double ll_out = has_outliers ? outlierComponent_ptr->outlier_likelihood(n) : -arma::datum::inf;
  const double log_w_non = has_outliers ? std::log(outlierComponent_ptr->non_outlier_weight) : 0.0;
  const double log_w_out = has_outliers ? std::log(outlierComponent_ptr->outlier_weight) : -arma::datum::inf;
  
  // Observed-data log-likelihood of the item under this view's own normalised 
  // weights (the cross-view coupling is not a property of a single view)
  const vec log_pi = log_weights - logSumExp(log_weights);
  double observed = logSumExp(log_pi + ll);
  if(has_outliers) {
    observed = logSumExp(vec({log_w_non + observed, log_w_out + ll_out}));
  }
  observed_likelihood_vec(n) = observed;
  
  if(fixed(n) == 0) {
    // Joint draw of the component and the outlier status. Component and 
    // outlier status are sampled together, so each conditional is exact.
    vec log_prob(has_outliers ? 2 * K : K);
    log_prob.subvec(0, K - 1) = log_component_weight + ll_sampling + log_w_non;
    if(has_outliers) {
      log_prob.subvec(K, 2 * K - 1) = log_component_weight + ll_out + log_w_out;
    }
    
    const double log_norm = logSumExp(log_prob);
    if(!std::isfinite(log_norm)) {
      Rcpp::stop("Non-finite allocation probabilities for item %d. Check for extreme "
                 "values in the data or degenerate component parameters.", (int) n + 1);
    }
    const vec prob = exp(log_prob - log_norm);
    
    const uword draw = sampleCategorical(prob);
    labels(n) = draw % K;
    outliers(n) = (draw >= K) ? 1 : 0;
    
    // Marginal probability of each component
    vec marginal = prob.subvec(0, K - 1);
    if(has_outliers) {
      marginal += prob.subvec(K, 2 * K - 1);
    }
    alloc.row(n) = marginal.t();
  }
  
  complete_likelihood_vec(n) = (outliers(n) == 1) ? ll_out : ll(labels(n));
}

void mixtureModel::sampleAllMissingValues() {
  const bool has_outliers = outlierComponent_ptr->active();
  density_ptr->labels = labels;
  for(uword n = 0; n < N; n++) {
    const arma::uvec& miss_idx = density_ptr->missing_indices(n);
    if(miss_idx.n_elem == 0) {
      continue;
    }
    if(has_outliers && outliers(n) == 1 && fixed(n) == 0) {
      const arma::vec sampled = outlierComponent_ptr->sampleMissingValues(n);
      for(uword i = 0; i < miss_idx.n_elem; i++) {
        density_ptr->X(n, miss_idx(i)) = sampled(i);
      }
    } else {
      density_ptr->sampleMissingForObservation(n);
    }
  }
}

void mixtureModel::initialiseMixture(
    const arma::vec& log_weights,
    const arma::mat& log_upweights
) {
  complete_likelihood = 0.0;
  observed_likelihood = 0.0;
  
  density_ptr->labels = labels;
  const vec log_pi = log_weights - logSumExp(log_weights);
  
  for (uword n = 0; n < N; n++) {
    const vec ll = density_ptr->itemLogLikelihood(n);
    observed_likelihood_vec(n) = logSumExp(log_pi + ll);
    complete_likelihood_vec(n) = ll(labels(n));
    
    if(fixed(n) == 0) {
      const vec log_prob = ll + log_weights + log_upweights.col(n);
      alloc.row(n) = exp(log_prob - logSumExp(log_prob)).t();
    }
  }
  observed_likelihood = accu(observed_likelihood_vec);
  complete_likelihood = accu(complete_likelihood_vec);
  { uvec occupied = unique(labels); K_occ = occupied.n_elem; }
}

void mixtureModel::setBeta(double beta_new) {
  if(!(beta_new >= 0.0 && beta_new <= 1.0)) {
    Rcpp::stop("The inverse temperature must lie in [0, 1].");
  }
  const bool has_missing = density_ptr->has_missing.n_elem > 0 && accu(density_ptr->has_missing) > 0;
  if(beta_new != 1.0 && (outlierComponent_ptr->active() || has_missing || mixture_type == 3)) {
    Rcpp::stop("Tempering supports complete data, no outlier component and "
               "'G', 'MVN' or 'C' densities only.");
  }
  beta = beta_new;
  density_ptr->beta = beta_new;
}

void mixtureModel::swapComponents(uword k, uword kprime) {
  density_ptr->swapComponents(k, kprime);
  alloc.swap_cols(k, kprime);
}


void mixtureModel::refreshLikelihoods(const arma::vec& log_weights) {
  const bool has_outliers = outlierComponent_ptr->active();
  const vec log_pi = log_weights - logSumExp(log_weights);
  for(uword n = 0; n < N; n++) {
    const vec ll = density_ptr->itemLogLikelihood(n);
    double observed = logSumExp(log_pi + ll);
    if(has_outliers) {
      observed = logSumExp(vec({std::log(outlierComponent_ptr->non_outlier_weight) + observed,
                                std::log(outlierComponent_ptr->outlier_weight) + outlierComponent_ptr->outlier_likelihood(n)}));
    }
    observed_likelihood_vec(n) = observed;
    complete_likelihood_vec(n) = (outliers(n) == 1) ? outlierComponent_ptr->outlier_likelihood(n) : ll(labels(n));
  }
  observed_likelihood = accu(observed_likelihood_vec);
  complete_likelihood = accu(complete_likelihood_vec);
}

bool mixtureModel::splitMergeMove(
    uword a, uword b,
    const arma::vec& log_weights,
    const arma::mat& log_upweights
) {
  const double beta_now = beta;
  const uword comp[2] = {a, b};
  const bool has_outliers = outlierComponent_ptr->active();
  const double log_w_non = has_outliers ? std::log(outlierComponent_ptr->non_outlier_weight) : 0.0;
  const double log_w_out = has_outliers ? std::log(outlierComponent_ptr->outlier_weight) : -arma::datum::inf;

  // The state of the move includes the imputed values of items with missing entries: a
  // non-outlier's complete-data density is part of its component's collapsed marginal, an
  // outlier's is the outlier law of the complete (observed + imputed) vector. The
  // sampler draws the imputations from exactly these conditionals, so this is the joint
  // they target. Without missing entries the two outlier densities coincide.
  auto outlier_density = [&](uword n) {
    if(density_ptr->missing_indices(n).n_elem == 0) {
      return outlierComponent_ptr->outlier_likelihood(n);
    }
    return outlierComponent_ptr->completeLogDensity(density_ptr->X.row(n).t());
  };

  // The free items of the two components and the statistics of the fixed ones
  // (an item flagged as an outlier does not enter the statistics of its component)
  collapsedStats base[2] = {density_ptr->emptyStats(), density_ptr->emptyStats()};
  std::vector<uword> items;
  for(uword n = 0; n < N; n++) {
    for(uword j = 0; j < 2; j++) {
      if(labels(n) == comp[j]) {
        if(fixed(n) == 1) {
          if(outliers(n) == 0) {
            density_ptr->addItemToStats(base[j], n);
          }
        } else {
          items.push_back(n);
        }
      }
    }
  }
  const uword m = items.size();
  if(m == 0) {
    return false;
  }

  // A uniformly random order (Fisher-Yates)
  for(uword i = m - 1; i > 0; i--) {
    const uword j = std::min<uword>((uword) std::floor(randu() * (double) (i + 1)), i);
    std::swap(items[i], items[j]);
  }

  // Walk the items in order, following either a random allocation (the proposal) or the 
  // current state, and return the sum of the log normalisers of the allocation 
  // probabilities. Each item chooses among (component, outlier flag): a non-outlier adds 
  // the collapsed predictive ratio of the component, an outlier adds the fixed outlier
  // density and leaves the component's statistics unchanged. A path entry is 2 * j + flag.
  auto walk = [&](bool follow_current, std::vector<uword>& path) {
    collapsedStats stats[2] = {base[0], base[1]};
    double lm[2] = {density_ptr->logMarginalLikelihood(stats[0], beta_now), 
                    density_ptr->logMarginalLikelihood(stats[1], beta_now)};
    double total = 0.0;
    path.assign(m, 0);
    for(uword t = 0; t < m; t++) {
      const uword n = items[t];
      collapsedStats cand[2] = {stats[0], stats[1]};
      double lm_cand[2];
      double score[4];
      for(uword j = 0; j < 2; j++) {
        density_ptr->addItemToStats(cand[j], n);
        lm_cand[j] = density_ptr->logMarginalLikelihood(cand[j], beta_now);
        const double base_score = log_weights(comp[j]) + log_upweights(comp[j], n);
        score[2 * j] = base_score + log_w_non + (lm_cand[j] - lm[j]);
        score[2 * j + 1] = has_outliers ? base_score + log_w_out + outlier_density(n)
                                        : -arma::datum::inf;
      }
      double mx = score[0];
      for(uword o = 1; o < 4; o++) {
        mx = std::max(mx, score[o]);
      }
      double acc = 0.0;
      for(uword o = 0; o < 4; o++) {
        acc += std::exp(score[o] - mx);
      }
      const double log_norm = mx + std::log(acc);
      if(!std::isfinite(log_norm)) {
        Rcpp::stop("Non-finite allocation probabilities in the split-merge move.");
      }
      total += log_norm;
      uword pick = 0;
      if(follow_current) {
        pick = 2 * ((labels(n) == comp[1]) ? 1 : 0) + ((outliers(n) == 1) ? 1 : 0);
      } else {
        const double u = randu();
        double cum = 0.0;
        pick = 3;
        for(uword o = 0; o < 4; o++) {
          cum += std::exp(score[o] - log_norm);
          if(u < cum) {
            pick = o;
            break;
          }
        }
        if(score[pick] == -arma::datum::inf) {
          // rounding at the upper end of the cumulative sum: take the last feasible option
          for(uword o = 4; o-- > 0;) {
            if(std::isfinite(score[o])) {
              pick = o;
              break;
            }
          }
        }
      }
      path[t] = pick;
      if(pick % 2 == 0) {
        stats[pick / 2] = cand[pick / 2];
        lm[pick / 2] = lm_cand[pick / 2];
      }
    }
    return total;
  };

  std::vector<uword> proposed_path, current_path;
  const double log_w_proposed = walk(false, proposed_path);
  const double log_w_current = walk(true, current_path);
  const double log_accept = log_w_proposed - log_w_current;
  const bool accept = (std::log(randu()) < log_accept);
  if(accept) {
    for(uword t = 0; t < m; t++) {
      labels(items[t]) = comp[proposed_path[t] / 2];
      outliers(items[t]) = proposed_path[t] % 2;
    }
    non_outliers = 1 - outliers;
  }
  return accept;
}
