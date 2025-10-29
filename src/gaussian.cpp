// gaussian.cpp
// =============================================================================
// included dependencies
# include "logLikelihoods.h"
# include "gaussian.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// gaussian class

gaussian::gaussian(arma::uword _K, arma::uvec _labels, arma::mat _X) : 
  density(_K, _labels, _X) 
{
  
  // Set the size of the objects to hold the component specific parameters
  mu.set_size(P, K);
  mu.zeros();
  
  std_devs.set_size(P, K);
  std_devs.zeros();
  
  // These will hold vertain matrix operations to avoid computational burden
  // This will hold the inverse of the squared standard devaitions
  precisions.set_size(P, K);
  precisions.zeros();
  
  log_std_devs.set_size(P, K);
  log_std_devs.zeros();
  
  // Mean vector and a diagonal covariance matrix and a component weight
  n_param = 2 * P + 1;
  
  // Default values for hyperparameters
  // Cluster hyperparameters for the Normal-inverse Wishart
  // Prior shrinkage
  kappa = 0.01;
  
  // Degrees of freedom
  nu = 3.0;
  
  // Empirical Bayesian hyperparameters for the mean and covariance
  empiricalBayesHyperparameters();
  
  identifyMissingValues();
  initializeMissingValues();
  
};


arma::vec gaussian::empiricalMean() {
  arma::vec mu_0(P);
  
  // Compute mean for each column using only finite values
  for(arma::uword p = 0; p < P; p++) {
    arma::vec col_data = X.col(p);
    arma::uvec finite_indices = arma::find_finite(col_data);
    
    if(finite_indices.n_elem > 0) {
      mu_0(p) = arma::mean(col_data.elem(finite_indices));
    } else {
      // If no finite values in this column, default to 0
      mu_0(p) = 0.0;
    }
  }
  
  return mu_0;
};

arma::mat gaussian::empiricalScaleVector() {
  double scale_entry = 0.0;
  arma::vec diag_entries(P);
  
  // Compute variances using only complete observations
  arma::uvec complete_obs;
  arma::mat X_complete;
  
  // Find rows with all finite values
  arma::uvec has_complete = arma::zeros<arma::uvec>(N);
  for(arma::uword n = 0; n < N; n++) {
    if(arma::all(arma::find_finite(X.row(n)))) {
      has_complete(n) = 1;
    }
  }
  complete_obs = arma::find(has_complete);
  
  if(complete_obs.n_elem > 1) {
    // Use complete observations for covariance
    X_complete = X.rows(complete_obs);
    arma::mat global_cov_loc = arma::cov(X_complete);
    
    // The entries of the diagonal of the empirical scale vector
    scale_entry = (arma::accu(global_cov_loc.diag()) / P) / std::pow(K, 2.0 / (double) P);
  } else {
    // Fallback: compute variance for each column separately using available data
    arma::vec col_vars(P);
    for(arma::uword p = 0; p < P; p++) {
      arma::vec col_data = X.col(p);
      arma::uvec finite_indices = arma::find_finite(col_data);
      
      if(finite_indices.n_elem > 1) {
        col_vars(p) = arma::var(col_data.elem(finite_indices));
      } else {
        col_vars(p) = 1.0; // Default variance
      }
    }
    scale_entry = arma::mean(col_vars) / std::pow(K, 2.0 / (double) P);
  }
  
  // Fill the vector
  diag_entries.fill(scale_entry);
  
  return diag_entries;
};

// arma::vec gaussian::empiricalMean() {
//   arma::vec mu_0;
//   arma::mat mean_mat;
//   mean_mat = arma::mean(X, 0).t();
//   mu_0 = mean_mat.col(0);
//   return mu_0;
// };
// 
// arma::mat gaussian::empiricalScaleVector() {
//   double scale_entry = 0.0;
//   arma::vec diag_entries(P);
//   arma::mat scale_param, global_cov_loc, Psi;
//   
//   
//   // Empirical Bayes for a diagonal covariance matrix
//   scale_param = X.each_row() - xi.t();
//   global_cov_loc = arma::cov(X);
//   
//   // The entries of the diagonal of the empirical scale matrix all have this 
//   // value
//   scale_entry = (arma::accu(global_cov_loc.diag()) / P) / std::pow(K, 2.0 / (double) P);
//   
//   // Fill the vector that corresponds to the diagonal entries of the scale matrix
//   diag_entries.fill(scale_entry);
// 
//   return diag_entries;
// };

void gaussian::empiricalBayesHyperparameters() {
  xi = empiricalMean();
  scale = empiricalScaleVector();
}

void gaussian::sampleStdDevPrior() {
  for(arma::uword k = 0; k < K; k++){
    for(uword p = 0; p < P; p++) {
      precisions(p, k) = randg(distr_param(0.5 * nu, 1.0 / (0.5 * scale(p))));
      std_devs(p, k) = 1.0 / precisions(p, k);
      log_std_devs(p, k) = std::log(std_devs(p, k));
    }
  }
};

void gaussian::sampleMuPrior() {
  for(arma::uword k = 0; k < K; k++){
    for(uword p = 0; p < P; p++) {
      mu(p, k) = randn() * (std_devs(p, k) / kappa) + xi(p);
    }
  }
};

void gaussian::sampleFromPriors() {
  sampleStdDevPrior();
  sampleMuPrior();
};

// The log likelihood of a item belonging to each cluster.
arma::vec gaussian::itemLogLikelihood(arma::uword n) {
  arma::vec ll(K);
  for(uword k = 0; k < K; k++) {
    ll(k) = logLikelihood(n, k);
  }
  return ll;
}

// // The log likelihood of a item belonging to a specific cluster.
// double gaussian::logLikelihood(arma::uword n, arma::uword k) {
//   arma::uvec obs_idx = observed_indices(n);
//   
//   if(obs_idx.n_elem == P) {
//     // Complete data - use existing fast calculation
//     arma::vec x_full = X.row(n).t();
//     arma::vec mu_k = mu.row(k).t();
//     arma::vec precision_k = precisions.row(k).t();
//     
//     // Diagonal covariance likelihood: sum of independent normal densities
//     double ll = 0.0;
//     for(uword p = 0; p < P; p++) {
//       double diff = x_full(p) - mu_k(p);
//       ll += -0.5 * (std::log(2.0 * M_PI) - std::log(precision_k(p)) + precision_k(p) * diff * diff);
//     }
//     return ll;
//     
//   } else if(obs_idx.n_elem > 0) {
//     // Missing data - sum over observed dimensions only
//     double ll = 0.0;
//     for(uword i = 0; i < obs_idx.n_elem; i++) {
//       uword p = obs_idx(i);
//       double x_val = X(n, p);
//       double mu_val = mu(k, p);
//       double precision_val = precisions(k, p);
//       
//       double diff = x_val - mu_val;
//       ll += -0.5 * (std::log(2.0 * M_PI) - std::log(precision_val) + precision_val * diff * diff);
//     }
//     return ll;
//   } else {
//     // All missing
//     return 0.0;
//   }
// }

// In gaussian.cpp - use direct indexing for diagonal covariance
double gaussian::logLikelihood(arma::uword n, arma::uword k) {
  arma::uvec obs_idx = observed_indices(n);
  
  if(obs_idx.n_elem == P) {
    // Complete data - use existing calculation
    arma::vec x_full = X.row(n).t();
    arma::vec mu_k = mu.col(k);
    arma::vec precision_k = precisions.col(k);
    
    double ll = 0.0;
    for(uword p = 0; p < P; p++) {
      double diff = x_full(p) - mu_k(p);
      ll += -0.5 * (std::log(2.0 * M_PI) - std::log(precision_k(p)) + precision_k(p) * diff * diff);
    }
    return ll;
  } else if(obs_idx.n_elem > 0) {
    // Missing data - loop over observed indices only
    double ll = 0.0;
    for(uword i = 0; i < obs_idx.n_elem; i++) {
      uword p = obs_idx(i);
      double x_val = X(n, p);           // Direct access
      double mu_val = mu(p, k);         // Direct access
      double precision_val = precisions(p, k);  // Direct access
      
      double diff = x_val - mu_val;
      ll += -0.5 * (std::log(2.0 * M_PI) - std::log(precision_val) + precision_val * diff * diff);
    }
    return ll;
  } else {
    return 0.0;
  }
}

void gaussian::sampleKthComponentParameters(
    uword k,
    umat members,
    uvec non_outliers
) {
  
  uword n_k = 0;
  double dist_from_prior = 0.0,
    kappa_n = 0.0, 
    nu_n = 0.0,
    scale_np = 0.0;
  
  uvec rel_inds;
  arma::vec mu_n(P), sample_mean(P), dist_from_mean;
  
  mat arma_cov(P, P), component_data, diff_from_mean;
  
  // Find the items relevant to sampling the parameters
  rel_inds = find((members.col(k) == 1) && (non_outliers == 1));
  
  // Find how many labels have the value
  n_k = rel_inds.n_elem;
  
  if(n_k > 0){
    
    // The vector that hold the distance of each observation from the mean of
    // the component data
    // dist_from_mean.reset();
    dist_from_mean.set_size(n_k);
    dist_from_mean.zeros();
    
    // component_data.reset();
    component_data.set_size(n_k, P);
    component_data.zeros();
    
    // diff_from_mean.reset();
    diff_from_mean.set_size(n_k, P);
    diff_from_mean.zeros();
    
    // Component data
    component_data = X.rows( rel_inds ) ;
    
    // Sample mean in the component data
    sample_mean = mean(component_data).t();
    
    mu_n = (xi * kappa + (double) n_k * sample_mean) / (kappa + (double) n_k);
    
    // Rcpp::Rcout << "\nIs this the issue?";
    diff_from_mean = component_data.each_row() - sample_mean.t();
    
    // Rcpp::Rcout << "\nEntering internal loop.";
    for(uword p = 0; p < P; p++) {
      kappa_n = kappa + (double) n_k;
      nu_n = nu + (double) n_k;
      
      // Rcpp::Rcout << "\nDist from mean.";
      dist_from_mean = arma::pow(diff_from_mean.col(p), 2.0);
      
      
      // Rcpp::Rcout << "\nDist from prior.";
      // Calculate the distance of the sample mean from the prior
      dist_from_prior = std::pow(sample_mean(p) - xi(p), 2.0);
      
      // Update the scale hyperparameter
      scale_np = (
        scale(p) 
        + accu(dist_from_mean) 
        + ((double) n_k * kappa / (kappa_n)) * dist_from_prior
      );
      
      // if( (0.5 * nu_n < 1e-6) || (1.0 / (0.5 * scale_np) < 1e-6) ) { 
      //   Rcpp::Rcout << "\nGaussian standarddeviation parameters below safe threshold.\n";
      //   Rcpp::Rcout << "\nnu: " << nu_n;
      //   Rcpp::Rcout << "\nscale_np: " << scale_np;
      //   Rcpp::Rcout << "\nReciprocal of scale_np: " << 1.0 / scale_np;
      // }
      
      // Sample the new precision
      precisions(p, k) = randg(distr_param(0.5 * nu_n, 1.0 / (0.5 * scale_np)));
      std_devs(p, k) = 1.0 / precisions(p, k);
      log_std_devs(p, k) = std::log(std_devs(p, k));
      
      // sample the new component mean in this measurement
      mu(p, k) = (randn() *  std_devs(p, k) / kappa_n) + mu_n(p);
    }
    // Rcpp::Rcout << "\nSampled mean:\n" << mu.col(k).t();
    // Rcpp::Rcout << "\nSampled std dev:\n" << std_devs.col(k).t();
    
    
    
  } else{
    // If no data in this component resample the parameters from the prio distn
    for(uword p = 0; p < P; p++) {
      
      precisions(p, k) = randg(distr_param(0.5 * nu, 1.0 / (0.5 * scale(p))));
      std_devs(p, k) = 1.0 / precisions(p, k);
      log_std_devs(p, k) = std::log(std_devs(p, k));
      
      mu(p, k) = randn() * (std_devs(p, k) / kappa) + xi(p);
    }
  }
};

void gaussian::initializeMissingValues() {
  // Same as MVN - use column means
  for(uword n = 0; n < N; n++) {
    if(missing_indices(n).n_elem > 0) {
      arma::uvec miss_idx = missing_indices(n);
      for(uword idx : miss_idx) {
        arma::vec col_data = X.col(idx);
        arma::uvec finite_indices = arma::find_finite(col_data);
        if(finite_indices.n_elem > 1) {
          double col_mean = arma::mean(col_data.elem(finite_indices));
          double col_sd = arma::stddev(col_data.elem(finite_indices));
          X(n, idx) = col_mean + arma::randn() * col_sd * 0.5;  // ADD JITTER
        } else if(finite_indices.n_elem == 1) {
          X(n, idx) = col_data(finite_indices(0)) + arma::randn() * 0.1;
        } else {
          X(n, idx) = arma::randn();
        }
      }
    }
  }
  X_t = X.t();
}

void gaussian::sampleMissingForObservation(arma::uword n) {
  if(missing_indices(n).n_elem > 0) {
    uword k = labels(n);
    arma::uvec miss_idx = missing_indices(n);
    
    // Independent sampling for diagonal covariance
    for(uword idx : miss_idx) {
      X(n, idx) = arma::randn() * std_devs(idx, k) + mu(idx, k);
    }
  }
}