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
  
};


arma::vec gaussian::empiricalMean() {
  arma::vec mu_0;
  arma::mat mean_mat;
  mean_mat = arma::mean(X, 0).t();
  mu_0 = mean_mat.col(0);
  return mu_0;
};

arma::mat gaussian::empiricalScaleVector() {
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

  return diag_entries;
};

void gaussian::empiricalBayesHyperparameters() {
  xi = zeros< vec >(P); // empiricalMean();
  scale = ones< vec >(P); // empiricalScaleVector();
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
arma::vec gaussian::itemLogLikelihood(arma::vec item) {
  arma::vec ll(K);
  ll.zeros();
  for(arma::uword k = 0; k < K; k++){
    ll(k) = logLikelihood(item, k);
  }
  return(ll);
};

// The log likelihood of a item belonging to a specific cluster.
double gaussian::logLikelihood(arma::vec item, arma::uword k) {
  
  double exponent = 0.0, ll = 0.0, dist_to_mean = 0.0;
  
  // The exponent part of the gaussian pdf
  for(uword p = 0; p < P; p++) {
    // ll += pNorm(item(p), mu(p, k), std_devs(p, k));
    dist_to_mean = std::pow(item(p) - mu(p, k), 2.0);
    exponent = dist_to_mean * precisions(p, k);

    // Normal log likelihood
    ll -= 0.5 *(log_std_devs(p, k) + exponent); 
  }
  ll -= 0.5 * (double) P * log(2.0 * M_PI);
  return(ll);
};

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

// void gaussian::sampleParameters(arma::umat members, arma::uvec non_outliers) {
// 
//   // for (arma::uword k = 0; k < K; k++) {
//   std::for_each(
//     std::execution::par,
//     K_inds.begin(),
//     K_inds.end(),
//     [&](uword k) {
//       sampleKthComponentParameters(k, members, non_outliers);
//     }
//   );
// };
