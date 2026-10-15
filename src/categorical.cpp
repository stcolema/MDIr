// categorical.cpp
// =============================================================================
// included dependencies
# include "logLikelihoods.h"
# include "categorical.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// categorical class

categorical::categorical(arma::uword _K, arma::uvec _labels, arma::mat _X) : 
  density(_K, _labels, _X) 
{
  n_cat.set_size(P);
  cat_prior_probability.set_size(P);
  category_probabilities.set_size(P);
  
  Y.zeros(N, P);
  
  identifyMissingValues();
  initializeMissingValues();
  initialiseParameters();
};

void categorical::initialiseParameters() {
  
  // Frequencies come from the observed entries only
  for(uword p = 0; p < P; p++) {
    const uvec obs = find(has_missing.col(p) == 0);
    uword max_cat = 0;
    for(uword i = 0; i < obs.n_elem; i++) {
      max_cat = std::max(max_cat, Y(obs(i), p));
    }
    n_cat(p) = max_cat + 1;
    
    category_probabilities(p).zeros(n_cat(p), K);
    
    vec freq(n_cat(p), arma::fill::zeros);
    for(uword i = 0; i < obs.n_elem; i++) {
      freq(Y(obs(i), p)) += 1.0;
    }
    freq /= std::max((double) obs.n_elem, 1.0);
    
    // Categories absent from the data (a gap in the coding) still need a
    // strictly positive concentration
    cat_prior_probability(p) = arma::clamp(freq, 1e-3, arma::datum::inf);
  } 
  
  n_param = sum(n_cat) * K;
}

Rcpp::List categorical::hyperparameterList() const {
  return Rcpp::List::create(
    Rcpp::Named("n_cat") = n_cat,
    Rcpp::Named("concentration") = cat_prior_probability
  );
}

void categorical::normaliseColumns(uword p) {
  for(uword k = 0; k < K; k++) {
    category_probabilities(p).col(k) /= accu(category_probabilities(p).col(k));
  }
}

void categorical::sampleFromPriors() {
  for(uword p = 0; p < P; p++) {
    for(uword k = 0; k < K; k++) {
      for(uword ii = 0; ii < n_cat(p); ii++) {
        category_probabilities(p)(ii, k) = rGamma(cat_prior_probability(p)(ii), 1.0);
      }
    }
    normaliseColumns(p);
  }
};

void categorical::sampleKthComponentParameters(
    uword k, 
    const umat& members, 
    const uvec& non_outliers
)  {
  const uvec relevant_indices = find((members.col(k) == 1) && (non_outliers == 1));
  
  for(uword p = 0; p < P; p++) {
    vec counts(n_cat(p), arma::fill::zeros);
    for(uword i = 0; i < relevant_indices.n_elem; i++) {
      counts(Y(relevant_indices(i), p)) += 1.0;
    }
    for(uword ii = 0; ii < n_cat(p); ii++) {
      category_probabilities(p)(ii, k) = rGamma(cat_prior_probability(p)(ii) + counts(ii), 1.0);
    }
    category_probabilities(p).col(k) /= accu(category_probabilities(p).col(k));
  }
}

arma::vec categorical::itemLogLikelihood(arma::uword n) {
  arma::vec ll(K);
  for(uword k = 0; k < K; k++) {
    ll(k) = logLikelihood(n, k);
  }
  return ll;
}

double categorical::logLikelihood(arma::uword n, arma::uword k) {
  const arma::uvec& obs_idx = observed_indices(n);
  double ll = 0.0;
  for(uword i = 0; i < obs_idx.n_elem; i++) {
    const uword p = obs_idx(i);
    // Probabilities are strictly positive up to underflow of a Gamma draw
    ll += std::log(std::max(category_probabilities(p)(Y(n, p), k), 
                            std::numeric_limits<double>::min()));
  }
  return ll;
}

void categorical::initializeMissingValues() {
  // Fill each missing entry with the modal observed category of its column
  for(uword p = 0; p < P; p++) {
    const uvec obs = find(has_missing.col(p) == 0);
    for(uword i = 0; i < obs.n_elem; i++) {
      const double value = X(obs(i), p);
      if(value < 0.0 || std::abs(value - std::round(value)) > 1e-8) {
        Rcpp::stop("Categorical data must be non-negative integers (column %d).", (int) p + 1);
      }
      Y(obs(i), p) = (uword) std::llround(value);
    }
    uword mode_val = 0;
    if(obs.n_elem > 0) {
      uword max_cat = 0;
      for(uword i = 0; i < obs.n_elem; i++) {
        max_cat = std::max(max_cat, Y(obs(i), p));
      }
      uvec counts(max_cat + 1, arma::fill::zeros);
      for(uword i = 0; i < obs.n_elem; i++) {
        counts(Y(obs(i), p))++;
      }
      mode_val = counts.index_max();
    }
    for(uword n = 0; n < N; n++) {
      if(has_missing(n, p) == 1) {
        Y(n, p) = mode_val;
        X(n, p) = (double) mode_val;
      }
    }
  }
}

void categorical::sampleMissingForObservation(arma::uword n) {
  const uword k = labels(n);
  const arma::uvec& miss_idx = missing_indices(n);
  for(uword i = 0; i < miss_idx.n_elem; i++) {
    const uword p = miss_idx(i);
    const uword sampled_category = sampleCategorical(category_probabilities(p).col(k));
    X(n, p) = (double) sampled_category;
    Y(n, p) = sampled_category;
  }
}

void categorical::swapComponents(uword k, uword kprime) {
  for(uword p = 0; p < P; p++) {
    category_probabilities(p).swap_cols(k, kprime);
  }
}

arma::vec categorical::parameters() const {
  arma::vec theta(n_param);
  uword start = 0;
  for(uword p = 0; p < P; p++) {
    const uword len = n_cat(p) * K;
    theta.subvec(start, start + len - 1) = vectorise(category_probabilities(p));
    start += len;
  }
  return theta;
}

void categorical::setParameters(const arma::vec& theta) {
  if(theta.n_elem != n_param) {
    Rcpp::stop("categorical: parameter vector has the wrong length.");
  }
  uword start = 0;
  for(uword p = 0; p < P; p++) {
    const uword len = n_cat(p) * K;
    category_probabilities(p) = reshape(theta.subvec(start, start + len - 1), n_cat(p), K);
    start += len;
  }
}

arma::vec categorical::simulate(arma::uword k) const {
  arma::vec x(P);
  for(uword p = 0; p < P; p++) {
    x(p) = (double) sampleCategorical(category_probabilities(p).col(k));
  }
  return x;
}

void categorical::replaceData(const arma::mat& X_new) {
  density::replaceData(X_new);
  Y.zeros(N, P);
  for(uword n = 0; n < N; n++) {
    const arma::uvec& obs_idx = observed_indices(n);
    for(uword i = 0; i < obs_idx.n_elem; i++) {
      const uword p = obs_idx(i);
      const double value = X(n, p);
      if(value < 0.0 || std::abs(value - std::round(value)) > 1e-8) {
        Rcpp::stop("Categorical data must be non-negative integers (column %d).", (int) p + 1);
      }
      const uword category = (uword) std::llround(value);
      if(category >= n_cat(p)) {
        Rcpp::stop("Column %d holds category %d, which the model was not fitted with (it has %d categories).",
                   (int) p + 1, (int) category, (int) n_cat(p));
      }
      Y(n, p) = category;
    }
  }
}
