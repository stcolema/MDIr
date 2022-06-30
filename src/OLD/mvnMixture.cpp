

# include "logLikelihoods.h"
# include "mvnMixture.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

mvnMixture::mvnMixture(
    arma::uword _K,
    arma::uvec _labels,
    // arma::vec _concentration,
    arma::mat _X
  ) : mixture(_K,
  _labels,
  // _concentration,
  _X) {

    // Mean vector and covariance matrix and a component weight
    n_param = (1 + P + P * (P + 1) * 0.5);

    // Default values for hyperparameters
    // Cluster hyperparameters for the Normal-inverse Wishart
    // Prior shrinkage
    kappa = 0.01;
    // Degrees of freedom
    nu = P + 2;

    // Mean
    arma::mat mean_mat = arma::mean(_X, 0).t();
    xi = mean_mat.col(0);

    // Empirical Bayes for a diagonal covariance matrix
    arma::mat scale_param = _X.each_row() - xi.t();
    arma::vec diag_entries(P);
    // double scale_entry = arma::accu(scale_param % scale_param, 0) / (N * std::pow(K, 1.0 / (double) P));

    arma::mat global_cov = arma::cov(X);
    double scale_entry = (arma::accu(global_cov.diag()) / P) / std::pow(K, 2.0 / (double) P);

    diag_entries.fill(scale_entry);
    scale = arma::diagmat( diag_entries );

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
  };


void mvnMixture::sampleCovPrior() {
  for(arma::uword k = 0; k < K; k++){
    cov.slice(k) = arma::iwishrnd(scale, nu);
    cov_inv.slice(k) = arma::inv_sympd(cov.slice(k));
    cov_log_det(k) = arma::log_det(cov.slice(k)).real();
  }
};

void mvnMixture::sampleMuPrior() {
  for(arma::uword k = 0; k < K; k++){
    mu.col(k) = arma::mvnrnd(xi, (1.0/kappa) * cov.slice(k), 1);
  }
};

void mvnMixture::sampleFromPriors() {
  sampleCovPrior();
  sampleMuPrior();
};

// Update the common matrix manipulations to avoid recalculating N times
void mvnMixture::matrixCombinations() {
  for(arma::uword k = 0; k < K; k++) {
    cov_inv.slice(k) = arma::inv_sympd(cov.slice(k));
    cov_log_det(k) = arma::log_det(cov.slice(k)).real();
  }
};

// The log likelihood of a item belonging to each cluster.
arma::vec mvnMixture::itemLogLikelihood(arma::vec item) {

  double exponent = 0.0;
  arma::vec ll(K), dist_to_mean(P);
  ll.zeros();
  dist_to_mean.zeros();

  for(arma::uword k = 0; k < K; k++){

    // The exponent part of the MVN pdf
    dist_to_mean = item - mu.col(k);
    exponent = arma::as_scalar(dist_to_mean.t() * cov_inv.slice(k) * dist_to_mean);

    // Normal log likelihood
    ll(k) = -0.5 *(cov_log_det(k) + exponent + (double) P * log(2.0 * M_PI));
  }
  return(ll);
};

// The log likelihood of a item belonging to a specific cluster.
double mvnMixture::logLikelihood(arma::vec item, arma::uword k) {
  
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

void mvnMixture::calcBIC(){

  // BIC = 2 * model_likelihood;

  BIC = 2 * complete_likelihood - n_param * K_occ * std::log(N);

  // for(arma::uword k = 0; k < K; k++) {
  //   BIC -= n_param_cluster * std::log(N_k(k) + 1);
  // }

};


void mvnMixture::sampleParameters() {

  arma::uword n_k = 0;
  arma::vec mu_n(P), sample_mean(P);
  arma::mat sample_cov(P, P), dist_from_prior(P, P), scale_n(P, P);

  for (arma::uword k = 0; k < K; k++) {

    // std::cout << "\nN_k (wrong): " << accu(labels == k);
    // std::cout << "\nN_k (should be right): " << N_k(k);


    // Find how many labels have the value
    n_k = N_k(k);
    
    // std::cout << "\n\nMembers:\n" << size(find(members.col(k)));
    // std::cout << "\n\nOutliers:\n" << size(find(non_outliers == 1));
    // std::cout << "\n\nNon-outlier members:\n" << size(find(members.col(k) && non_outliers == 1));
    
    if(n_k > 0){

      // std::cout << "\n\nSampling new value";
      
      // Component data
      arma::mat component_data = X.rows( arma::find(members.col(k) && (non_outliers == 1)) );

      // std::cout << "\n\nComponent data subsetted.\n" << size(component_data);
      
      // Sample mean in the component data
      sample_mean = arma::mean(component_data).t();

      // std::cout << "\n\nSampled mean acquired:\n" << sample_mean;
      
      sample_cov = calcSampleCov(component_data, sample_mean, n_k, P);

      // std::cout << "\n\nSampled cov acquired:\n" << sample_cov;
      
      // Calculate the distance of the sample mean from the prior
      dist_from_prior = (sample_mean - xi) * (sample_mean - xi).t();
      
      // std::cout << "\n\nDistance from prior calculated:\n" << dist_from_prior;
      
      // Update the scale hyperparameter
      scale_n = scale + sample_cov + ((kappa * n_k) / (double) (kappa + n_k)) * dist_from_prior;
      
      // std::cout << "\n\nDistance from prior calculated:\n" << dist_from_prior;
      // 
      // std::cout << "\n\nCovariacne:\n" << cov.slice(k);
      // std::cout << "\n\nNu: " << nu;
      // std::cout << "\nN_k: " << n_k;
      
      // Sample a new covariance matrix
      cov.slice(k) = arma::iwishrnd(scale_n, nu + n_k);
      
      // std::cout << "\n\nCovariance sampled";
      
      // The weighted average of the prior mean and sample mean
      mu_n = (kappa * xi + n_k * sample_mean) / (double)(kappa + n_k);
      
      // std::cout << "\nMu_n calculated";
      
      // Sample a new mean vector
      mu.col(k) = arma::mvnrnd(mu_n, (1.0 / (double) (kappa + n_k)) * cov.slice(k), 1);
      
      // std::cout << "\nMean sampled";
      // std::cout << "\n\nDistance from prior calculated:\n" << dist_from_prior;
      
    } else{

      // If no members in the component, draw from the prior distribution
      cov.slice(k) = arma::iwishrnd(scale, nu);
      mu.col(k) = arma::mvnrnd(xi, (1.0 / (double) kappa) * cov.slice(k), 1);

    }
    
    // Save the inverse and log determinant of the new covariance matrices
    cov_inv.slice(k) = arma::inv_sympd(cov.slice(k));
    cov_log_det(k) = arma::log_det(cov.slice(k)).real();
    
  }
};

// double mvnMixture::marginalLikelihood(uvec indices) {
//   
//   mat component_data = X.rows(indices);
//   
//   uword n_k = length(indices);
//   arma::vec mu_n(P), sample_mean(P);
//   arma::mat sample_cov(P, P), dist_from_prior(P, P), scale_n(P, P);
//   
//   // Sample mean in the component data
//   sample_mean = arma::mean(component_data).t();
//   sample_cov = calcSampleCov(component_data, sample_mean, n_k, P);
//       
//   // Calculate the distance of the sample mean from the prior
//   dist_from_prior = (sample_mean - xi) * (sample_mean - xi).t();
//   
//   // Update the scale hyperparameter
//   scale_n = scale + sample_cov + ((kappa * n_k) / (double) (kappa + n_k)) * dist_from_prior;
//       
//   
// }

double mvnMixture::posteriorPredictive(arma::vec x, arma::uvec indices) {

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
