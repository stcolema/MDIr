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

// The log likelihood of a item belonging to each cluster.
arma::vec mvn::itemLogLikelihood(arma::vec item) {
  arma::vec ll(K), dist_to_mean(P);
  ll.zeros();
  dist_to_mean.zeros();
  
  for(arma::uword k = 0; k < K; k++){
    // Normal log likelihood
    ll(k) = logLikelihood(item, k);
  }
  return(ll);
};

// The log likelihood of a item belonging to a specific cluster.
double mvn::logLikelihood(arma::vec item, arma::uword k) {
  
  double exponent = 0.0, ll = 0.0;
  arma::vec dist_to_mean(P);
  dist_to_mean.zeros();
  
  // The exponent part of the MVN pdf
  dist_to_mean = item - mu.col(k);
  exponent = arma::as_scalar(dist_to_mean.t() * cov_inv.slice(k) * dist_to_mean);
  
  // Normal log likelihood
  ll = -0.5 *(cov_log_det(k) + exponent + (double) P * log(2.0 * M_PI));
  
  return(ll);
};
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
