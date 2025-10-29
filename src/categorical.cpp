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
  // category_probabilities.zeros();
  
  Y = conv_to<umat>::from(X);
  
  identifyMissingValues();
  initializeMissingValues();
  
  // Initialise some of the more awkward parameters
  initialiseParameters();

};


void categorical::initialiseParameters() {
  uvec Y_p(N), categories;
  mat call_prob_entry;
  
  // Rcpp::Rcout << "\nInitialise parameters in categorical densities.";
  
  for(uword p = 0; p < P; p++) {
    
    // Rcpp::Rcout << "\nAccess pth column.";
    
    Y_p = Y.col(p);
    
    // Rcpp::Rcout << "\nFind the number of categories in each measurement/feature.";
    
    // Find the number of categories in the column
    categories = unique(Y_p); 
    n_cat(p) = categories.n_elem;

    // Rcpp::Rcout << "\nDefine the entry in the class probabilities field.";
        
    // Create a matrix of 0's. This is a placeholder for the probability for 
    // each measurement within each cluster.
    call_prob_entry.set_size(n_cat(p), K);
    call_prob_entry.zeros();
    
    // Rcpp::Rcout << "\nSet the entry to this.";
    category_probabilities(p) = call_prob_entry;
    
    // Set the prior probability of being in any category empirically
    // Rcpp::Rcout << "\nSet the prior probability empricially.";
    cat_prior_probability(p).set_size(n_cat(p));
    for(uword ii = 0; ii < n_cat(p); ii++) {
      cat_prior_probability(p)(ii) = ((double) accu(Y_p == ii)) / (double) N;
    }
    
    // Reset the entry to empty.
    call_prob_entry.reset();
  } 
  
  n_param = sum(n_cat) * K;
  
}

void categorical::sampleFromPriors() {
  
  // Rcpp::Rcout << "\nSample from categorical prior distribution.";
  
  for(uword p = 0; p < P; p++) {
    for(uword ii = 0; ii < n_cat(p); ii++) {

      category_probabilities(p).row(ii) = rGamma(K, cat_prior_probability(p)(ii), 1.0).t();
        //   arma::randg(
        //   K,
        //   arma::distr_param(cat_prior_probability(p)(ii), 1.0)
        // ).t();
      // }
    }
    
    for(uword k = 0; k < K; k++) {
      category_probabilities(p).col(k) *= 1.0 / accu(category_probabilities(p).col(k));
    }
  }
};

void categorical::sampleKthComponentParameters(
    uword k, 
    umat members, 
    uvec non_outliers
)  {
  uvec relevant_indices;
  umat component_data;
  uword cat_count = 0;
  double concentration_n = 0.0;
  
  // Find the items relevant to sampling the parameters
  relevant_indices = find((members.col(k) == 1) && (non_outliers == 1));
  
  component_data = Y.rows(relevant_indices);
  for(uword p = 0; p < P; p++) {
    
    for(uword ii = 0; ii < n_cat(p); ii++) {
      cat_count = accu(component_data.col(p) == ii);
      
      concentration_n = cat_prior_probability(p)(ii) + cat_count;
      
      category_probabilities(p)(ii, k) = rGamma(concentration_n, 1.0);
      // arma::randg(
      //   K, 
      //   arma::distr_param(concentration_n, 1.0)
      // ).t();
      
    }
    category_probabilities(p).col(k) *= 1.0 / accu(category_probabilities(p).col(k));
  }
}

// void categorical::sampleParameters(arma::umat members, arma::uvec non_outliers) {
//   uvec relevant_indices;
//   umat component_data;
//   uword cat_count = 0;
//   double concentration_n = 0.0;
//   
//   std::for_each(
//     std::execution::par,
//     K_inds.begin(),
//     K_inds.end(),
//     [&](uword k) {
//       sampleKthComponentParameters(k, members, non_outliers);
//     }
//   );
//   
//   // for(uword k = 0; k < K; k++) {
//   //   
//   //   // Find the items relevant to sampling the parameters
//   //   relevant_indices = find((members.col(k) == 1) && (non_outliers == 1));
//   //   
//   //   component_data = Y.rows(relevant_indices);
//   //   for(uword p = 0; p < P; p++) {
//   //     
//   //     for(uword ii = 0; ii < n_cat(p); ii++) {
//   //       cat_count = accu(component_data.col(p) == ii);
//   //       
//   //       concentration_n = cat_prior_probability(p)(ii) + cat_count;
//   //       
//   //       category_probabilities(p).row(ii) = arma::randg(
//   //         K, 
//   //         arma::distr_param(concentration_n, 1.0)
//   //       ).t();
//   //     }
//   //   }
//   // }
// };

// double categorical::logLikelihood(arma::vec item, arma::uword k) {
//   
//   double ll = 0.0;
//   uword x_p = 0;
//   
//   for(uword p = 0; p < P; p++) {
//     x_p = item(p);
//     ll += std::log(category_probabilities(p)(x_p, k));
//   }
//   
//   return ll;
// };
// 
// arma::vec categorical::itemLogLikelihood(arma::vec item) {
//   vec ll(K);
//   
//   for(uword k = 0; k < K; k++) {
//     ll(k) = logLikelihood(item, k);
//   }
//   return ll;
// };

arma::vec categorical::itemLogLikelihood(arma::uword n) {
  arma::vec ll(K);
  for(uword k = 0; k < K; k++) {
    ll(k) = logLikelihood(n, k);
  }
  return ll;
}

double categorical::logLikelihood(arma::uword n, arma::uword k) {
  arma::uvec obs_idx = observed_indices(n);
  
  if(obs_idx.n_elem == P) {
    // Complete data - use all features
    double ll = 0.0;
    for(uword p = 0; p < P; p++) {
      uword category = Y(n, p);  // Integer category value
      
      // Safety check for valid category
      if(category < category_probabilities(p).n_rows) {
        double prob = category_probabilities(p)(category, k);
        // Add small epsilon to avoid log(0)
        ll += std::log(prob + 1e-10);
      } else {
        // Invalid category - assign very low probability
        ll += std::log(1e-10);
      }
    }
    return ll;
    
  } else if(obs_idx.n_elem > 0) {
    // Missing data - use observed features only
    double ll = 0.0;
    for(uword i = 0; i < obs_idx.n_elem; i++) {
      uword p = obs_idx(i);
      uword category = Y(n, p);
      
      if(category < category_probabilities(p).n_rows) {
        double prob = category_probabilities(p)(category, k);
        ll += std::log(prob + 1e-10);
      } else {
        ll += std::log(1e-10);
      }
    }
    return ll;
  } else {
    // All missing
    return 0.0;
  }
}


void categorical::initializeMissingValues() {
  for(uword n = 0; n < N; n++) {
    if(missing_indices(n).n_elem > 0) {
      arma::uvec miss_idx = missing_indices(n);
      for(uword idx : miss_idx) {
        arma::vec col_data = X.col(idx);
        arma::uvec finite_indices = arma::find_finite(col_data);
        if(finite_indices.n_elem > 0) {
          // Find mode (most frequent category)
          arma::vec observed_values = col_data.elem(finite_indices);
          arma::vec unique_vals = arma::unique(observed_values);
          uword mode_val = 0;
          uword max_count = 0;
          for(uword i = 0; i < unique_vals.n_elem; i++) {
            uword count = arma::sum(observed_values == unique_vals(i));
            if(count > max_count) {
              max_count = count;
              mode_val = static_cast<uword>(unique_vals(i));
            }
          }
          X(n, idx) = static_cast<double>(mode_val);
        } else {
          X(n, idx) = 0.0;
        }
      }
    }
  }
  
  // Update integer matrix Y
  for(uword n = 0; n < N; n++) {
    for(uword p = 0; p < P; p++) {
      Y(n, p) = static_cast<uword>(X(n, p));
    }
  }
  X_t = X.t();
}

void categorical::sampleMissingForObservation(arma::uword n) {
  if(missing_indices(n).n_elem > 0) {
    uword k = labels(n);
    arma::uvec miss_idx = missing_indices(n);
    
    for(uword idx : miss_idx) {
      arma::vec probs = category_probabilities(idx).col(k);
      double u = arma::randu();
      arma::vec cumprobs = arma::cumsum(probs);
      uword sampled_category = 0;
      for(uword cat = 0; cat < cumprobs.n_elem; cat++) {
        if(u <= cumprobs(cat)) {
          sampled_category = cat;
          break;
        }
      }
      X(n, idx) = static_cast<double>(sampled_category);
      Y(n, idx) = sampled_category;
    }
  }
}