// mvn.cpp
// =============================================================================
// included dependencies
# include "logLikelihoods.h"
# include "mvn.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// mvn class

mvn::mvn(arma::uword _K, arma::uvec _labels, arma::mat _X) : 
  density(_K, _labels, _X) 
{
  
  // Set the size of the objects to hold the component specific parameters
  mu.set_size(P, K);
  mu.zeros();
  
  cov.set_size(P, P, K);
  cov.zeros();
  
  // These will hold vertain matrix operations to avoid computational burden
  // The log determinant of each cluster covariance
  cov_log_det = arma::zeros<arma::vec>(K);
  
  // Inverse of the cluster covariance
  cov_inv.set_size(P, P, K);
  cov_inv.zeros();
  
  // Mean vector and covariance matrix and a component weight
  n_param = P * (1 + (P + 1) * 0.5);
  
  // Default values for hyperparameters
  // Cluster hyperparameters for the Normal-inverse Wishart
  // Prior shrinkage
  kappa = 0.01;
  // Degrees of freedom
  nu = P + 2;
  
  // Empirical Bayesian hyperparameters for the mean and covariance
  empiricalBayesHyperparameters();
  
  // Identify and initialize the missing values
  identifyMissingValues();
  initializeMissingValues();
};


arma::vec mvn::empiricalMean() {
  arma::vec mu_0;
  arma::mat mean_mat;
  mean_mat = arma::mean(X, 0).t();
  mu_0 = mean_mat.col(0);
  return mu_0;
};

arma::mat mvn::empiricalScaleMatrix() {
  double scale_entry = 0.0;
  arma::vec diag_entries(P);
  arma::mat scale_param, global_cov_loc, Psi;
  
  
  // Empirical Bayes for a diagonal covariance matrix
  scale_param = X.each_row() - xi.t();
  global_cov_loc = arma::cov(X);
  
  // The entries of the diagonal of the empirical scale matrix all have this 
  // value
  scale_entry = (arma::accu(global_cov_loc.diag()) / P) / std::pow(K, 2.0 / (double) P);
  
  // Fill the vector that corresponds to the diagonal entries of the scale matrix
  diag_entries.fill(scale_entry);
  
  // The empirical scale matrix
  Psi = arma::diagmat( diag_entries );
  return Psi;
};

void mvn::empiricalBayesHyperparameters() {
  xi = empiricalMean();
  scale = empiricalScaleMatrix();
}

void mvn::sampleCovPrior() {
  // Rcpp::Rcout << "\nScale:\n" << scale;
  // Rcpp::Rcout << "\nDF: " << nu << "\n";
  for(arma::uword k = 0; k < K; k++){
    cov.slice(k) = arma::iwishrnd(scale, nu);
    cov_inv.slice(k) = arma::inv_sympd(cov.slice(k));
    cov_log_det(k) = arma::log_det_sympd(cov.slice(k));
  }
  // Rcpp::Rcout << "\nCovariances sampled from prior.\n";
};

void mvn::sampleMuPrior() {
  for(arma::uword k = 0; k < K; k++){
    mu.col(k) = arma::mvnrnd(xi, (1.0/kappa) * cov.slice(k), 1);
  }
};

void mvn::sampleFromPriors() {
  sampleCovPrior();
  sampleMuPrior();
  matrixCombinations();
};

// Update the common matrix manipulations to avoid recalculating N times
void mvn::matrixCombinations() {
  for(arma::uword k = 0; k < K; k++) {
    cov_inv.slice(k) = arma::inv_sympd(cov.slice(k));
    cov_log_det(k) = arma::log_det_sympd(cov.slice(k));
  }
};

// // Modified likelihood function
// double mvn::logLikelihood(arma::uword n, arma::uword k) {
//   arma::uvec obs_idx = observed_indices(n);
//   
//   if(obs_idx.n_elem == P) {
//     // Complete data - fast path
//     return pNorm(X.row(n).t(), mu.col(k), cov.slice(k), true);
//   } else {
//     // Missing data - marginal likelihood on observed dimensions
//     arma::vec x_obs = X.row(n).elem(obs_idx).t();
//     arma::vec mu_obs = mu.col(k).elem(obs_idx);
//     arma::mat cov_obs = cov.slice(k).submat(obs_idx, obs_idx);
//     return pNorm(x_obs, mu_obs, cov_obs, true);
//   }
// }

// // The log likelihood of a item belonging to a specific cluster.
// double mvn::logLikelihood(arma::vec item, arma::uword k) {
//   
//   double exponent = 0.0, ll = 0.0;
//   arma::vec dist_to_mean(P);
//   dist_to_mean.zeros();
//   
//   // The exponent part of the MVN pdf
//   dist_to_mean = item - mu.col(k);
//   exponent = arma::as_scalar(dist_to_mean.t() * cov_inv.slice(k) * dist_to_mean);
//   
//   // Normal log likelihood
//   ll = -0.5 *(cov_log_det(k) + exponent + (double) P * log(2.0 * M_PI));
//   
//   return(ll);
// };

void mvn::sampleParameters(arma::umat members, arma::uvec non_outliers) {
  
  // for(uword k = 0; k < K; k++) {
  std::for_each(
    std::execution::par,
    K_inds.begin(),
    K_inds.end(),
    [&](uword k) {
      sampleKthComponentParameters(k, members, non_outliers);
    }
  );
  
  matrixCombinations();
};


void mvn::sampleKthComponentParameters(
    uword k, 
    umat members, 
    uvec non_outliers
  ) {
  
  arma::uword n_k = 0;
  uvec rel_inds;
  arma::vec mu_n(P), sample_mean(P);
  arma::mat sample_cov(P, P), dist_from_prior(P, P), scale_n(P, P), component_data;
  mat arma_cov(P, P);
  
  // Find the items relevant to sampling the parameters
  rel_inds = find((members.col(k) == 1) && (non_outliers == 1));
  
  // Find how many labels have the value
  n_k = rel_inds.n_elem;
  
  if(n_k > 0){
    
    // Component data
    component_data = X.rows( rel_inds ) ;
  
    // Sample mean in the component data
    sample_mean = sampleMean(component_data);
    
    // Sample covariance times its degree of freedom
    sample_cov = calcSampleCov(component_data, sample_mean, n_k, P);
    // arma_cov = (n_k - 1) * arma::cov(component_data);
    
    // Calculate the distance of the sample mean from the prior
    dist_from_prior = (sample_mean - xi) * (sample_mean - xi).t();
    
    // Update the scale hyperparameter
    scale_n = scale + sample_cov + ((kappa * (double) n_k) / (kappa + (double) n_k)) * dist_from_prior;
    
    // Sample a new covariance matrix
    cov.slice(k) = iwishrnd(scale_n, nu + (double) n_k);
    
    // The weighted average of the prior mean and sample mean
    mu_n = (kappa * xi + (double) n_k * sample_mean) / (kappa + (double) n_k);
    
    // Sample a new mean vector
    mu.col(k) = mvnrnd(mu_n, (1.0 / (kappa + (double) n_k)) * cov.slice(k));
    
  } else{
    
    // If no members in the component, draw from the prior distribution
    cov.slice(k) = iwishrnd(scale, nu);
    mu.col(k) = mvnrnd(xi, (1.0 / (double) kappa) * cov.slice(k));
    
  }

  // Save the inverse and log determinant of the new covariance matrices
  cov_inv.slice(k) = inv_sympd(cov.slice(k));
  cov_log_det(k) = log_det_sympd(cov.slice(k));

};

double mvn::posteriorPredictive(arma::vec x, arma::uvec indices) {
  
  mat component_data = X.rows(indices);
  
  uword n_k = indices.n_rows;
  double nu_n_rel = nu + n_k - P + 1, kappa_n = kappa + n_k;
  arma::vec mu_n(P), sample_mean(P);
  arma::mat sample_cov(P, P), dist_from_prior(P, P), scale_n(P, P);
  
  // Sample mean in the component data
  sample_mean = arma::mean(component_data).t();
  
  mu_n = (kappa * xi + n_k * sample_mean) / (double)(kappa + n_k);
  
  sample_cov = calcSampleCov(component_data, sample_mean, n_k, P);
  
  // Calculate the distance of the sample mean from the prior
  dist_from_prior = (sample_mean - xi) * (sample_mean - xi).t();
  
  // Update the scale hyperparameter
  scale_n = scale + sample_cov + ((kappa * n_k) / (double) (kappa + n_k)) * dist_from_prior;
  
  return mvtLogLikelihood(x, mu_n, scale_n / (kappa_n * nu_n_rel), nu_n_rel);
};



void mvn::initializeMissingValues() {
  for(uword n = 0; n < N; n++) {
    if(missing_indices(n).n_elem > 0) {
      arma::uvec miss_idx = missing_indices(n);
      for(uword idx : miss_idx) {
        arma::vec col_data = X.col(idx);
        arma::uvec finite_indices = arma::find_finite(col_data);
        if(finite_indices.n_elem > 0) {
          X(n, idx) = arma::mean(col_data.elem(finite_indices));
        } else {
          X(n, idx) = arma::randn() * 0.1;
        }
      }
    }
  }
  X_t = X.t();
}

void mvn::sampleMissingForObservation(arma::uword n) {
  if(missing_indices(n).n_elem > 0) {
    uword k = labels(n);
    arma::uvec miss_idx = missing_indices(n);
    arma::uvec obs_idx = observed_indices(n);
    
    if(obs_idx.n_elem > 0) {
      // CORRECT: Get full vectors first, then extract elements
      arma::vec x_n_full = X.row(n).t();              // Full observation n as column vector
      arma::vec mu_k_full = mu.col(k);            // Full component k parameters as column vector
      
      arma::vec x_obs = x_n_full.elem(obs_idx);       // Extract observed elements
      arma::vec mu_obs = mu_k_full.elem(obs_idx);     // Extract observed parameters
      arma::vec mu_miss = mu_k_full.elem(miss_idx);   // Extract missing parameters
      
      arma::mat cov_miss = cov.slice(k).submat(miss_idx, miss_idx);
      arma::mat cov_obs = cov.slice(k).submat(obs_idx, obs_idx);
      arma::mat cov_cross = cov.slice(k).submat(miss_idx, obs_idx);
      
      arma::vec conditional_mean = mu_miss + cov_cross * arma::solve(cov_obs, x_obs - mu_obs);
      arma::mat conditional_cov = cov_miss - cov_cross * arma::solve(cov_obs, cov_cross.t());
      conditional_cov = 0.5 * (conditional_cov + conditional_cov.t());
      
      // Numerical safety
      arma::vec eigval = arma::eig_sym(conditional_cov);
      if(eigval.min() < 1e-8) {
        conditional_cov += 1e-6 * arma::eye(conditional_cov.n_rows, conditional_cov.n_cols);
      }
      
      arma::vec sampled = arma::mvnrnd(conditional_mean, conditional_cov);
      
      // Update missing values
      for(uword i = 0; i < miss_idx.n_elem; i++) {
        X(n, miss_idx(i)) = sampled(i);
      }
    } else {
      // All missing
      arma::vec mu_k_full = mu.col(k);
      arma::vec mu_miss = mu_k_full.elem(miss_idx);
      arma::mat cov_miss = cov.slice(k).submat(miss_idx, miss_idx);
      arma::vec sampled = arma::mvnrnd(mu_miss, cov_miss);
      for(uword i = 0; i < miss_idx.n_elem; i++) {
        X(n, miss_idx(i)) = sampled(i);
      }
    }
  }
}

double mvn::logLikelihood(arma::uword n, arma::uword k) {
  arma::uvec obs_idx = observed_indices(n);
  
  if(obs_idx.n_elem == P) {
    // Complete data
    return pNorm(X.row(n).t(), mu.col(k), cov.slice(k), true);
  } else if(obs_idx.n_elem > 0) {
    // CORRECT: Get full vectors first, then extract elements
    arma::vec x_n_full = X.row(n).t();
    arma::vec mu_k_full = mu.col(k);
    
    arma::vec x_obs = x_n_full.elem(obs_idx);
    arma::vec mu_obs = mu_k_full.elem(obs_idx);
    arma::mat cov_obs = cov.slice(k).submat(obs_idx, obs_idx);
    
    return pNorm(x_obs, mu_obs, cov_obs, true);
  } else {
    return 0.0;
  }
}

arma::vec mvn::itemLogLikelihood(arma::uword n) {
  arma::vec ll(K);
  for(uword k = 0; k < K; k++) {
    ll(k) = logLikelihood(n, k);
  }
  return ll;
}