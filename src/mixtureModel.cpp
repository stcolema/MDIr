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

// BIC currently ignores outlier parameters
void mixtureModel::calcBIC() {
  BIC = 2 * complete_likelihood - (n_param + 1) * K_occ * std::log((double) N);
}

void mixtureModel::updateOutlierWeights() {
  non_outliers = 1 - outliers;
  outlierComponent_ptr->updateWeights(non_outliers, outliers);
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
    log_prob.subvec(0, K - 1) = log_component_weight + ll + log_w_non;
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

void mixtureModel::swapComponents(uword k, uword kprime) {
  density_ptr->swapComponents(k, kprime);
  alloc.swap_cols(k, kprime);
}
