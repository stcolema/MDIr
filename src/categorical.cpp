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

double categorical::logLikelihood(arma::vec item, arma::uword k) {
  
  double ll = 0.0;
  uword x_p = 0;
  
  for(uword p = 0; p < P; p++) {
    x_p = item(p);
    ll += std::log(category_probabilities(p)(x_p, k));
  }
  
  return ll;
};

arma::vec categorical::itemLogLikelihood(arma::vec item) {
  vec ll(K);
  
  for(uword k = 0; k < K; k++) {
    ll(k) = logLikelihood(item, k);
  }
  return ll;
};
