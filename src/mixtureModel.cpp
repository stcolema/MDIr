// mixtureModel.cpp
// =============================================================================
// included dependencies
# include <RcppParallel.h>
# include <RcppArmadillo.h>
# include <execution>
# include "mixtureModel.h"

using namespace arma ;

// =============================================================================
// mixtureModel class

// Parametrised class
mixtureModel::mixtureModel(
  arma::uword _mixture_type,
  arma::uword _outlier_type,
  arma::uword _K,
  arma::uvec _labels,
  arma::uvec _fixed,
  arma::mat _X) {
  
  mixture_type = _mixture_type;
  outlier_type = _outlier_type;
  K = _K;
  labels = _labels;
  X = _X;
  X_t = X.t();
  
  // Dimensions
  N = X.n_rows;
  P = X.n_cols;
  
  N_inds = linspace < uvec >(0, N - 1, N);
  
  // Class populations
  N_k = zeros<uvec>(K);
  
  // Log likelihood (individual and model)
  ll = zeros<vec>(K);
  complete_likelihood_vec = zeros< vec >(N);
  observed_likelihood_vec = zeros< vec >(N);
  likelihood = zeros<vec>(N);
  
  // Class members
  members.set_size(N, K);
  members.zeros();
  
  // Allocation probability (not sensible in unsupervised)
  alloc.set_size(N, K);
  alloc.zeros();
  
  // Pass the indicator vector for being fixed to the ``fixed`` object.
  fixed = _fixed;
  N_fixed = accu(fixed);
  fixed_ind = find(fixed == 1);
  unfixed_ind = find(fixed == 0);
  
  // Set the known label allocations to 1
  for(uword n = 0; n < N; n++) {
    if(fixed(n) == 1) {
      alloc(n, labels(n)) = 1.0;
    }
  }
  
  // Used a few times
  vec_of_ones = ones<uvec>(N);
  
  // Outlier vectors; not used really here
  non_outliers = ones<uvec>(N);
  outliers = zeros<uvec>(N);
  
  // Rcpp::Rcout << "\nInitialise density within mixture.\n";
  
  // Initialise the density
  initialiseDensity(mixture_type);
  
  // Rcpp::Rcout << "\nAccess density pointer.\n";
  
  n_param = density_ptr->n_param;
  
  
  // Rcpp::Rcout << "\nInitialise outlier component within mixture.\n";
  
  // Initialise the outlier component (default is empty)
  initialiseOutlierComponent(outlier_type);
  
  // Share missing patterns with outlier component after both are initialized
  outlierComponent_ptr->setMissingPatterns(
      density_ptr->missing_indices,
      density_ptr->observed_indices, 
      density_ptr->has_missing
  );
  
  // Initialize missing values in outlier component
  outlierComponent_ptr->initializeMissingValues();
  
  acceptance_count.set_size(0);
  hypers.set_size(0);
  if(mixture_type == 3) {
    acceptance_count.set_size(3 * K);
    acceptance_count.zeros();
    hypers.set_size(3 * K);
    hypers.zeros();
  }
  density_ptr->acceptance_count = acceptance_count;
  density_ptr->hypers = hypers;
};

// void mixtureModel::updateItemAllocation(
//     uword n, 
//     vec log_weights, 
//     vec log_upweigths
//   ) {
//   double u = 0.0;
//   uvec uniqueK;
//   vec comp_prob(K), ll(K);
//   
//   // The mixture-specific log likelihood for each observation in each class
//   ll = itemLogLikelihood(X_t.col(n));
//   
//   // Update with weights
//   comp_prob = ll + log_weights + log_upweigths;
//   
//   // Record the likelihood - this is used to calculate the observed likelihood
//   observed_likelihood_vec(n) = accu(comp_prob);
//   
//   if(fixed(n) == 0) {
//     // Handle overflow problems and then normalise to convert to probabilities
//     comp_prob = exp(comp_prob - max(comp_prob));
//     comp_prob = comp_prob / sum(comp_prob);
//     
//     // Save the allocation probabilities
//     alloc.row(n) = comp_prob.t();
//     
//     // Prediction and update
//     u = randu<double>( );
//     
//     labels(n) = sum(u > cumsum(comp_prob));
//     outliers(n) = sampleOutlier(ll(labels(n)), outlier_likelihood(n));
//   }
//   
//   // Update the complete likelihood based on the new labelling
//   // complete_likelihood += ll(labels(n));
//   complete_likelihood_vec(n) = ll(labels(n));
//   
// };
  
// void mixtureModel::updateAllocation(
//   arma::vec log_weights, 
//   arma::mat log_upweigths
// ) {
//   uvec uniqueK;
//   vec comp_prob(K);
//   
//   complete_likelihood = 0.0;
//   observed_likelihood = 0.0;
//   
//   // Update the outlier weights - if no outlier component is being modelled, 
//   // this does nothing
//   updateOutlierWeights();
//   
//   // for (auto& n : unfixed_ind) {
//   std::for_each(std::execution::par, N_inds.begin(), N_inds.end(), [&] (uword n) {
//     updateItemAllocation(n, log_weights, log_upweigths.col(n));
//   }
//   );
//   
//   observed_likelihood = accu(observed_likelihood_vec);
//   complete_likelihood = accu(complete_likelihood_vec);
//   
//   // Number of occupied components (used in BIC calculation)
//   uniqueK = unique(labels);
//   K_occ = uniqueK.n_elem;
// };

void mixtureModel::updateAllocation(arma::vec log_weights, arma::mat log_upweigths) {
  uvec uniqueK;
  vec comp_prob(K);
  
  complete_likelihood = 0.0;
  observed_likelihood = 0.0;
  
  // 1. Update outlier weights
  updateOutlierWeights();
  
  
  // 2. Sample missing values
  sampleAllMissingValues();
  
  // 3. Update cluster assignments using current missing value imputations
  std::for_each(std::execution::par, N_inds.begin(), N_inds.end(), [&] (uword n) {
    updateItemAllocation(n, log_weights, log_upweigths.col(n));
  });
  
  // 4. Aggregate likelihoods (will be recalculated next iteration)
  observed_likelihood = accu(observed_likelihood_vec);
  complete_likelihood = accu(complete_likelihood_vec);
  
  // 5. Update component counts
  uniqueK = unique(labels);
  K_occ = uniqueK.n_elem;
}  

void mixtureModel::updateItemAllocation(uword n, vec log_weights, vec log_upweigths) {
  double u = 0.0;
  uvec uniqueK;
  vec comp_prob(K), ll(K);
  
  // Use new signature for likelihood functions
  ll = density_ptr->itemLogLikelihood(n);
  
  // Update with weights
  comp_prob = ll + log_weights + log_upweigths;
  
  // Record the likelihood - this is used to calculate the observed likelihood
  observed_likelihood_vec(n) = accu(comp_prob);
  
  if(fixed(n) == 0) {
    // Handle overflow problems and then normalise to convert to probabilities
    comp_prob = exp(comp_prob - max(comp_prob));
    comp_prob = comp_prob / sum(comp_prob);
    
    // Save the allocation probabilities
    alloc.row(n) = comp_prob.t();
    
    // Prediction and update
    u = randu<double>( );
    
    labels(n) = sum(u > cumsum(comp_prob));
    outliers(n) = sampleOutlier(ll(labels(n)), outlier_likelihood(n));
  }
  
  // Update the complete likelihood based on the new labelling
  // complete_likelihood += ll(labels(n));
  complete_likelihood_vec(n) = ll(labels(n));
}

void mixtureModel::sampleAllMissingValues() {
  for(uword n = 0; n < N; n++) {
    if(density_ptr->missing_indices(n).n_elem > 0) {
      if(outliers(n) == 1 && fixed(n) == 0) {
        // Sample from outlier component
        outlierComponent_ptr->sampleMissingForObservation(n);
      } else {
        // Sample from mixture component
        density_ptr->sampleMissingForObservation(n);
      }
    }
  }
  // Single coordinated update of transpose
  X_t = X.t();
  density_ptr->X_t = X_t;
  outlierComponent_ptr->X_t = X_t;
}

void mixtureModel::initialiseDensity(arma::uword type) {
  
  densityFactory my_factory;
  
  // Convert the unsigned integer into a mixture type object
  densityFactory::densityType val = static_cast<densityFactory::densityType>(type);
  
  // Create a smart pointer to the correct type of model
  density_ptr = my_factory.createDensity(val, K, labels, X);
};
  
// Sample from the prior distribution of the density
void mixtureModel::sampleFromPriors() {
  density_ptr->sampleFromPriors();
};

void mixtureModel::sampleParameters() {
  density_ptr->sampleParameters(members, non_outliers);
};

// BIC currently ignores outlier parametes
void mixtureModel::calcBIC() {
  BIC = 2 * complete_likelihood - (n_param + 1) * K_occ * std::log(N);
}
// arma::vec mixtureModel::itemLogLikelihood(arma::vec x) {
//   return density_ptr->itemLogLikelihood(x);
// };

arma::vec mixtureModel::itemLogLikelihood(arma::uword n) {
  return density_ptr->itemLogLikelihood(n);
};

// double mixtureModel::logLikelihood(arma::vec x, arma::uword k) {
//   return density_ptr->logLikelihood(x, k);
// };

double mixtureModel::logLikelihood(arma::uword n, arma::uword k) {
  return density_ptr->logLikelihood(n, k);
};


void mixtureModel::initialiseOutlierComponent(uword type) {
  
  // Rcpp::Rcout << "\nIn function to initialise outlier component.\n";
  
  outlierComponentFactory my_factory;

  // Rcpp::Rcout << "\nCast type.\n";
  
  // Convert the unsigned integer into a mixture type object
  outlierComponentFactory::outlierType val = static_cast<outlierComponentFactory::outlierType>(type);

  // Rcpp::Rcout << "\nCreate pointer.\n";
  
  // Create a smart pointer to the correct type of model
  outlierComponent_ptr = my_factory.createOutlierComponent(val, fixed, X);
  
  // Rcpp::Rcout << "\nAccess likelihood.\n";
  outlier_likelihood = outlierComponent_ptr->outlier_likelihood;
  
  // Outlier vectors
  // Rcpp::Rcout << "Access outliers.\n";
  
  outliers = outlierComponent_ptr->outliers;
  
  // Rcpp::Rcout << "Update non-outliers.\n";
  non_outliers = 1 - outliers;
  
};


void mixtureModel::initialiseMixture(
    arma::vec log_weights,
    arma::mat log_upweigths
) {
  
  uvec uniqueK;
  vec comp_prob(K);
  
  complete_likelihood = 0.0;
  observed_likelihood = 0.0;
  
  // for (auto& n : unfixed_ind) {
  for (uword n = 0; n < N; n++) {
    
    // The mixture-specific log likelihood for each observation in each class
    // ll = itemLogLikelihood(X_t.col(n));
    ll = itemLogLikelihood(n);
    
    // Update with weights
    comp_prob = ll + log_weights + log_upweigths.col(n);
    
    // Record the likelihood - this is used to calculate the observed likelihood
    // likelihood(n) = accu(comp_prob);
    observed_likelihood += accu(comp_prob);
    
    if(fixed(n) == 0) {
      // Handle overflow problems and then normalise to convert to probabilities
      comp_prob = exp(comp_prob - max(comp_prob));
      comp_prob = comp_prob / sum(comp_prob);
      
      // Save the allocation probabilities
      alloc.row(n) = comp_prob.t();
    }
    
    // Update the complete likelihood based on the new labelling
    complete_likelihood += ll(labels(n));
    
  }
  
  // Number of occupied components (used in BIC calculation)
  uniqueK = unique(labels);
  K_occ = uniqueK.n_elem;
  
}

arma::uword mixtureModel::sampleOutlier(double non_outlier_likelihood_n,
                                        double outlier_likelihood_n) {
  return outlierComponent_ptr->sampleOutlier(non_outlier_likelihood_n, 
                                      outlier_likelihood_n);
};

void mixtureModel::updateOutlierWeights() {
  non_outliers = 1 - outliers;
  outlierComponent_ptr->updateWeights(non_outliers, outliers);
}

// void mixtureModel::sampleMissingData() {
//   umat data_mask = ones(N, P);
//   mat X_inferred = X;
//   field<umat> missing_idx(N);
//   
//   for(uword n = 0; n < N; n++) {
//     missing_idx.row(n) = find(data_mask.row(n) == 1);
//     
//   // std::for_each(std::execution::par, N_inds.begin(), N_inds.end(), [&] (uword n) {
//     X_inferred(n, missing_idx) = sampleMissingDataPoint(
//       X.row(n).t(), 
//       labels(n), 
//       missing_idx.row(n).t()
//     ).t();
//   // }
//   // );
//   }
// };
