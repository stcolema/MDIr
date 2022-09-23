// gp.cpp
// =============================================================================
// included dependencies
# include "logLikelihoods.h"
# include "gp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// mvn class

gp::gp(arma::uword _K, arma::uvec _labels, arma::mat _X) : 
  density(_K, _labels, _X) 
{
  
  zero_vec.set_size(P);
  zero_vec.zeros();
  
  // Hyperparameters
  amplitude.set_size(K);
  amplitude.ones();
  
  length.set_size(K);
  length.ones();
    
  noise.set_size(K);
  noise.ones();
 
  // Set the size of the objects to hold the component specific parameters
  mu.set_size(P, K);
  mu.zeros();
  
  kernel_sub_block.set_size(P, P, K);
  kernel_sub_block.zeros();
  
  time_diff_mat.set_size(P, P);
  time_diff_mat.zeros();
  
  for(int ii = 0; ii < P; ii++) {
    for(int jj = ii + 1; jj < P; jj++) {
      time_diff_mat(ii, jj) = - 0.5 * std::pow((double) (jj - ii), 2.0);
      time_diff_mat(jj, ii) = time_diff_mat(ii, jj);
    }
  }
  
  I_p = eye(P, P);
  
  noise_acceptance_count.set_size(K);
  noise_acceptance_count.zeros();
  
  length_acceptance_count.set_size(K);
  length_acceptance_count.zeros();
  
  amplitude_acceptance_count.set_size(K);
  amplitude_acceptance_count.zeros();
  
  // Those weird GP objects that can be massive and vary from component to 
  // component, the size being a function of N_k.
  repeated_time_indices.set_size(K);
  repeated_mean_vector.set_size(K);
  covariance_matrix.set_size(K);
  flattened_component_data.set_size(K);
  
  t_inds = arma::linspace< uvec >(0, P - 1, P);
  time_difference_mat.set_size(P, P);
  time_difference_mat.zeros();
  for(uword ii = 0; ii < (P - 1); ii++) {
    for(uword jj = (ii + 1); jj < P; jj++) {
      time_difference_mat(ii, jj) = jj - ii;
      time_difference_mat(jj, ii) = jj - ii;
    }
  }
  time_difference_mat = arma::pow(time_difference_mat, 2.0);
  
  // These will hold vertain matrix operations to avoid computational burden
  // The log determinant of each cluster covariance
  cov_log_det = arma::zeros<arma::vec>(K);
  
  // Inverse of the cluster covariance
  inverse_covariance.set_size(K);
  
  // Mean vector and covariance matrix and a component weight
  // n_param = P * (1 + (P + 1) * 0.5);
  
  // Empirical Bayesian hyperparameters for the mean and covariance
  // empiricalBayesHyperparameters();
  
  hypers.set_size(3 * K);
  hypers.zeros();
  
  acceptance_count.set_size(3 * K);
  acceptance_count.zeros();
  
};

double gp::noisePriorLogDensity(double x, bool logNorm) {
  double y = 0.0;
  if(logNorm) {
    y = pNorm(log(x), 0.0, noise_prior_std_dev);
  } else {
    y = pHalfCauchy(x, 0.0, 5.0, true);
  }
  return y;
}

double gp::ampltiduePriorLogDensity(double x, bool logNorm) {
  double y = 0.0;
  if(logNorm) {
    y = pNorm(log(x), 0.0, hyper_prior_std_dev);
  } else {
    y = pHalfCauchy(x, 0.0, 5.0, true);
  }
  return y;
}

double gp::lengthPriorLogDensity(double x, bool logNorm) {
  double y = 0.0;
  if(logNorm) {
    y = pNorm(log(x), 0.0, hyper_prior_std_dev);
  } else {
    y = pHalfCauchy(x, 0.0, 5.0, true);
  }
  return y;
}

double gp::sampleLengthPriorDistribution(bool logNorm, double threshold) {
  double x = 0.0;
  if(logNorm) {
    x = std::exp(randn< double >() * hyper_prior_std_dev);
  } else {
    x = rHalfCauchy(0.0, 5.0);
  }
  if(x < threshold) {
    x = sampleLengthPriorDistribution(logNorm, threshold);
  }
  return x;
};

double gp::sampleAmplitudePriorDistribution(bool logNorm, double threshold) {
  double x = 0.0;
  if(logNorm) {
    x = std::exp(randn< double >() * hyper_prior_std_dev);
  } else {
    x = rHalfCauchy(0.0, 5.0);
  }
  if(x < threshold) {
    x = sampleAmplitudePriorDistribution(logNorm, threshold);
  }
  return x;
};

double gp::sampleNoisePriorDistribution(bool logNorm, double threshold) {
  double x = 0.0;
  if(logNorm) {
    x = std::exp(randn< double >() * noise_prior_std_dev);
  } else {
    x = rHalfCauchy(0.0, 5.0);
  }
  // x = std::exp(randn< double >());
  if(x < threshold) {
    x = sampleNoisePriorDistribution(logNorm, threshold);
  }
  return x;
};

void gp::sampleKthComponentHyperParameterPrior(uword k, bool logNorm) {
  amplitude(k) = sampleAmplitudePriorDistribution(logNorm);
  length(k) = sampleLengthPriorDistribution(logNorm);
  noise(k) = sampleNoisePriorDistribution(logNorm);
};

void gp::sampleHyperParameterPriors() {
  for(uword k = 0; k < K; k++) {
    sampleKthComponentHyperParameterPrior(k, logNormPriorUsed);
  }
  hypers.subvec(0, K - 1) = amplitude;
  hypers.subvec(K, 2 * K - 1) = length;
  hypers.subvec(2 * K, 3 * K - 1) = noise;
};

// arma::mat gp::calculateCovarianceKernel(arma::uvec t_inds) {
//   return cov_kernel_ptr->calculateCovarianceMatrix(t_inds);
// };

// std::unique_ptr<kernel> gp::initialiseKernel(uword kernel_type) {
// 
//   kernelFactory my_factory;
//   
//   // Convert the unsigned integer into a mixture type object
//   kernelFactory::kernelType val = static_cast<kernelFactory::kernelType>(kernel_type);
//   
//   // Create a smart pointer to the correct type of model
//   std::unique_ptr<kernel> kernel_ptr = my_factory.createKernel(val);
// 
//   return kernel_ptr;
// }

// === Prior-related functions =================================================
void gp::sampleMuPrior() {
  for(arma::uword k = 0; k < K; k++){
    mu.col(k) = arma::mvnrnd(zero_vec, kernel_sub_block.slice(k));
  }
};

void gp::sampleFromPriors() {
  sampleHyperParameterPriors();
  calculateKernelSubBlock();
  sampleMuPrior();
};

// === Covariance function =====================================================

mat gp::calculateKthComponentKernelSubBlock(double amplitude,
                                            double length,
                                            double kernel_subblock_threshold) {
  mat sub_block(P, P);
  sub_block.zeros();
  
  sub_block = std::log(amplitude) + (1.0 / length) * time_diff_mat;
  sub_block = exp(sub_block);
  
  // for(uword ii = 0; ii < P; ii++) {
  //   sub_block(ii, ii) = amplitude;
  //   for(uword jj = ii + 1; jj < P; jj++) {
  //     sub_block(ii, jj) = squaredExponentialFunction(
  //       amplitude,
  //       length,
  //       ii,
  //       jj
  //     );
  // 
  //     if(sub_block(ii, jj) < kernel_subblock_threshold) {
  //       sub_block(ii, jj) = 0.0;
  //       sub_block(jj, ii) = 0.0;
  //       // break;
  //     }
  //     sub_block(jj, ii) = sub_block(ii, jj);
  //   }
  // }
  
  return sub_block;
};

void gp::calculateKernelSubBlock() {
  kernel_sub_block.ones();
  for(uword k = 0; k < K; k++) {
    kernel_sub_block.slice(k) = calculateKthComponentKernelSubBlock(
      amplitude(k),
      length(k)
    );
  }
};

mat gp::constructCovarianceMatrix(uword n_k, mat kernel_sub_block) {
  mat covariance_matrix(n_k * P, n_k * P);
  covariance_matrix.zeros();
  covariance_matrix = repmat(kernel_sub_block, n_k, n_k);
  return covariance_matrix;
};

mat gp::smallerInversion(uword n_k, double noise, mat kernel_sub_block) {
  mat Q = I_p + (((double) n_k) / noise) * kernel_sub_block;
  return inv_sympd(Q);
}

mat gp::firstCovProduct(uword n_k, double noise, mat kernel_sub_block) {
  mat B(P, P), Z(P, P), output(P, P);
  output.zeros();
  Z.zeros();
  B.zeros();
  
  Z = smallerInversion(n_k, noise, kernel_sub_block);
  B = I_p - Z;
  
  output = (1.0 / noise) * (kernel_sub_block - (kernel_sub_block * B));
  return output;
};


mat gp::invertComponentCovariance(uword n_k, double noise, mat kernel_sub_block) {
  mat J(n_k, n_k), Q_k(P, P), Z_k(P, P), I_NkP(n_k * P, n_k * P), out(n_k * P, n_k * P);
  J.ones();
  Q_k.zeros();
  Z_k.zeros();

  I_NkP = eye(n_k * P, n_k * P);
  Q_k = I_p + ( (double)  n_k / noise) * kernel_sub_block;
  // Z_k = smallerInversion(n_k, noise, kernel_sub_block);
  Z_k = inv_sympd(Q_k);
  out = (1.0 / noise) * I_NkP - (1.0 / ((double) n_k * noise)) * kron(J, I_p - Z_k);
  return out;
};


mat gp::covCheck(mat C, bool checkSymmetry, bool checkStability, double threshold) {
  bool not_symmetric = false, not_invertible = false, not_sympd = false;
  vec eigval(P);
  mat small_identity = 1e-7 * I_p;
  // C = roundMatrix(C, threshold);
  // not_sympd = ! C.is_sympd();
  // if(not_sympd) {
  //   Rcpp::Rcout << "\nNot symmetric positive definite.\n";
  // }
  // We can have that the covariance matrix becomes asymetric; this appears to 
  // be a floating point error, so we hardcode that the matrix is symmetric #
  // based on the uuper right traingle of the calculated covariance martix
  if(checkSymmetry) {
    
    mat u_cov = trimatu(C,  1);  // omit the main diagonal
    mat l_cov = trimatl(C, -1).t();  // omit the main diagonal
    
    not_symmetric = ! C.is_symmetric();
    // bool not_symmetric_my_check = approx_equal(u_cov, l_cov, "reldiff", 0.1);
    if(not_symmetric) {
      // Rcpp::Rcout << "\nNot symmetric. Reconstructing from upper right triangular matrix.\n";
      mat new_cov(P, P); // u_cov = trimatu(C, 1);
      new_cov = u_cov + u_cov.t();
      new_cov.diag() = C.diag();
      C = new_cov;
    }
  }
  // If our covariance matrix is poorly behaved (i.e. non-invertible), add a 
  // small constant to the diagonal entries
  if(checkStability) {
    eigval = eig_sym( C );
    not_invertible = min(eigval) < 1e-8;
    if(not_invertible) {
      // Rcpp::Rcout << "\nNot numerical stable for inversion. Add constant to diagonal.\n";
      small_identity *= 1e-7;
      C += small_identity;
    }
  }
  return C;
};

// mat gp::calculateCovTilde(uword n_k, double noise, mat cov_mat) {
//   mat first_product(P, P), final_product(P, P), cov_tilde(P, P);
//   first_product = firstCovProduct(n_k, noise, cov_mat);
//   final_product = ((double) n_k) * (first_product * cov_mat);
//   cov_tilde = cov_mat - final_product;
//   return cov_tilde;
// }

// === Mean function posterior==================================================

vec gp::posteriorMeanParameter(
    mat data, 
    mat first_product
) {
  uword n = data.n_rows;
  vec mu_tilde(P), sample_mean(P); 
  sample_mean = sampleMean(data);
  mu_tilde = (double) n * first_product * sample_mean;
  
  return mu_tilde;
};


// mat gp::posteriorCovarianceParameter(
//     mat covariance_matrix,
//     mat inverse_covariance_matrix) {
//   mat cov_tilde(P, P);
//   cov_tilde.zeros();
//   
//   cov_tilde = covariance_matrix.submat(P_inds, P_inds)
//     - (covariance_matrix.rows(P_inds)
//     * inverse_covariance_matrix
//     * covariance_matrix.cols(P_inds));
// 
//   return cov_tilde;
// };

vec gp::sampleMeanFunction(vec mu_tilde, mat cov_tilde) {
  vec eigval;
  mat eigvec, eigval_mat(P, P), stochasticity = mvnrnd(zeros<vec>(P), eye(P, P));
  eigval_mat.zeros();
  
  // Find negative eigen values and set to zero as these are not stable
  eig_sym( eigval, eigvec, cov_tilde );
  eigval.elem(find(eigval < 0.0) ).zeros();
  eigval_mat.diag() = arma::pow(eigval, 0.5);
  
  return mu_tilde + eigvec * eigval_mat * stochasticity;
};

void gp::sampleMeanPosterior(uword k, uword n_k, mat data) {
  bool sampleHypers = false;
  vec mu_tilde(P), sample_mean(P);
  mat
    cov_tilde(P, P), 
    covariance_matrix(n_k * P, n_k * P),
    inverse_covariance(n_k * P, n_k * P),
    rel_cov_mat(P, P),
    first_product(P, P),
    final_product(P, P);
  
  sample_mean = sampleMean(data);
  
  // Objects related to the covariance function
  rel_cov_mat = kernel_sub_block.slice(k); // covariance_matrix.rows(P_inds);

  // The product of the covariance matrix and the inverse as used in sampling 
  // parameters.
  first_product = firstCovProduct(n_k, noise(k), rel_cov_mat);
  final_product = ((double) n_k) * (first_product * rel_cov_mat);
  
  // Mean and covariance hyperparameter
  mu_tilde = ((double) n_k) * first_product * sample_mean;
  cov_tilde = rel_cov_mat - final_product;
  
  // Check that the covariance hyperparameter is numerically stable, add some 
  // small value to the diagonal if necessary
  // cov_tilde = covCheck(cov_tilde, false, true, matrix_precision);
  
  mu.col(k) = sampleMeanFunction(mu_tilde, cov_tilde);
  
  sampleHypers = (
    (samplingCount < 100 && (samplingCount % sampleHypersFrequencyBefore100) == 0) ||
    (samplingCount > 100 && samplingCount < 1000 && (samplingCount % sampleHypersFrequencyBefore1000) == 0) ||
    (samplingCount > 1000 && (samplingCount % sampleHypersFrequencyAfter1000) == 0)
  );
  
  if(sampleHypers) {
    sampleHyperParametersKthComponent(
      k,
      n_k,
      mu_tilde,
      sample_mean,
      cov_tilde
    );
    sampleNoise(k, n_k, data);
  }
};

void gp::sampleKthComponentParameters(uword k, umat members, uvec non_outliers) {
  
  // Find the items relevant to sampling the parameters
  uvec rel_inds = find((members.col(k) == 1) && (non_outliers == 1));
  mat component_data;
  
  // Find how many labels have the value
  uword n_k = rel_inds.n_elem;
  if(n_k > 0){
    // Component data
    component_data.set_size(n_k);
    component_data = X.rows( rel_inds ) ;
    
    // Sample parameters
    sampleMeanPosterior(k, n_k, component_data);
  } else {
    // Sample from the prior
    mu.col(k) = arma::mvnrnd(zero_vec, kernel_sub_block.slice(k));
    sampleKthComponentHyperParameterPrior(k);
  }
};

void gp::sampleParameters(arma::umat members, arma::uvec non_outliers) {
  arma::uword n_k = 0;
  uvec rel_inds;
  calculateKernelSubBlock();
  
  // for(uword k = 0; k < K; k++) {
  std::for_each(
    std::execution::par,
    K_inds.begin(),
    K_inds.end(),
    [&](uword k) {
      sampleKthComponentParameters(k, members, non_outliers);
    }
  );
  
  samplingCount++;
  
  hypers.subvec(0, K - 1) = amplitude;
  hypers.subvec(K, 2 * K - 1) = length;
  hypers.subvec(2 * K, 3 * K - 1) = noise;
  
  acceptance_count.subvec(0, K - 1) = amplitude_acceptance_count;
  acceptance_count.subvec(K, 2 * K - 1) = length_acceptance_count;
  acceptance_count.subvec(2 * K, 3 * K - 1) = noise_acceptance_count;
  
  // Rcpp::Rcout << "\n\n\nPARAMETERS";
  // Rcpp::Rcout << "\nMembership:\n" << sum(members);
  // Rcpp::Rcout << "\nNoise:\n" << noise.t(); 
  // Rcpp::Rcout << "\nNoise acceptance count:\n" << noise_acceptance_count.t();
  // 
  // Rcpp::Rcout << "\n\nLength:\n" << length.t();
  // Rcpp::Rcout << "\nLength acceptance count:\n" << length_acceptance_count.t();
  // 
  // Rcpp::Rcout << "\n\nAmplitude:\n" << amplitude.t();
  // Rcpp::Rcout << "\nAmplitude acceptance count:\n" << amplitude_acceptance_count.t();
  // 
  // Rcpp::Rcout << "\n\nMean:\n" << mu;
  
};

// === Hyper-parameters ========================================================

void gp::receiveHyperParametersProposalWindows(vec proposal_windows) {
  amplitude_proposal_window = proposal_windows[0];
  length_proposal_window = proposal_windows[1];
  noise_proposal_window = proposal_windows[2];
}

// Need to sample the hyperparameter and recalculate the mu tilde / cov tilde
double gp::hyperParameterLogKernel(
    double hyper, 
    vec mu_k, 
    vec mu_tilde, 
    mat cov_tilde, 
    bool logNorm
  ) {
  double score = 0.0;
  // Likelihood contribution
  score = pNorm(mu_k, mu_tilde, cov_tilde, false);
  // Prior contribution
  if(logNorm) {
    score += pNorm(log(hyper), 0.0, hyper_prior_std_dev);
  } else {
    score += pHalfCauchy(hyper, 0.0, 5.0);
  }
  return score;
};

void gp::sampleLength(
    uword k, 
    uword n_k, 
    vec mu_tilde, 
    vec sample_mean, 
    mat cov_tilde,
    double threshold
) {
  bool accept = false;
  double 
    acceptance_prob = 0.0, 
    new_score = 0.0, 
    old_score = 0.0,
    new_length = 0.0;
  vec new_mu_tilde(P);
  mat 
    new_sub_block(P, P), 
    new_cov_mat(n_k * P, n_k * P),
    new_inv_cov_mat(n_k * P, n_k * P), 
    new_cov_tilde(P, P),
    first_product(P, P),
    first_product_repeated(P, n_k * P),
    final_product(P, P);
  
  new_length = proposeNewNonNegativeValue(
    length(k), 
    length_proposal_window,
    use_log_norm_proposal,
    threshold
  );
  
  new_sub_block = calculateKthComponentKernelSubBlock(amplitude(k), new_length);
  
  // The product of the covariance matrix and the inverse as used in sampling 
  // parameters.
  first_product = firstCovProduct(n_k, noise(k), new_sub_block);
  final_product = (double) n_k * (first_product.cols(P_inds) * new_sub_block);
  
  // new_mu_tilde = first_product_repeated * component_data;
  new_mu_tilde = (double) n_k * first_product * sample_mean;
  new_cov_tilde = new_sub_block - final_product;
  
  // new_cov_tilde = covCheck(new_cov_tilde, false, true, matrix_precision);
  // 
  // if(rcond(new_cov_tilde) < threshold) {
  //   return;
  // }
  
  new_score = hyperParameterLogKernel(
    new_length, 
    mu.col(k), 
    new_mu_tilde, 
    new_cov_tilde,
    logNormPriorUsed
  );
  
  old_score = hyperParameterLogKernel(
    length(k), 
    mu.col(k), 
    mu_tilde, 
    cov_tilde,
    logNormPriorUsed
  );
  
  acceptance_prob =  std::min(1.0, std::exp(new_score - old_score));
  accept = metropolisAcceptanceStep(acceptance_prob);
  if(accept) {
    length(k) = new_length;
    length_acceptance_count(k)++;
  }
};

void gp::sampleAmplitude(
    uword k, 
    uword n_k, 
    vec mu_tilde, 
    vec sample_mean, 
    mat cov_tilde,
    double threshold
  ) {
  bool accept = false;
  double 
    acceptance_prob = 0.0, 
    new_score = 0.0, 
    old_score = 0.0,
    new_amplitude = 0.0;
  vec new_mu_tilde(P);
  mat 
    new_sub_block(P, P), 
    new_cov_mat(n_k * P, n_k * P),
    new_inv_cov_mat(n_k * P, n_k * P), 
    new_cov_tilde(P, P),
    first_product(P, P),
    first_product_repeated(P, n_k * P),
    final_product(P, P);
  
  new_amplitude = proposeNewNonNegativeValue(
    amplitude(k), 
    amplitude_proposal_window,
    use_log_norm_proposal
  );
    // std::exp(std::log(amplitude(k) + randn() * amplitude_proposal_window));
  if(new_amplitude < threshold) {
    return;
  }
  
  new_sub_block = calculateKthComponentKernelSubBlock(new_amplitude, length(k));
  first_product = firstCovProduct(n_k, noise(k), new_sub_block);
  final_product = (double) n_k * (first_product.cols(P_inds) * new_sub_block);

  new_mu_tilde = (double) n_k * first_product * sample_mean;
  new_cov_tilde = new_sub_block - final_product;
  
  // new_cov_tilde = covCheck(new_cov_tilde, false, true, matrix_precision);
  // 
  // if(rcond(new_cov_tilde) < threshold) {
  //   return;
  // }
  
  new_score = hyperParameterLogKernel(
    new_amplitude, 
    mu.col(k), 
    new_mu_tilde, 
    new_cov_tilde,
    logNormPriorUsed
  );
  
  old_score = hyperParameterLogKernel(
    amplitude(k), 
    mu.col(k), 
    mu_tilde, 
    cov_tilde,
    logNormPriorUsed
  );

  acceptance_prob =  std::min(1.0, std::exp(new_score - old_score));
  accept = metropolisAcceptanceStep(acceptance_prob);
  if(accept) {
    amplitude(k) = new_amplitude;
    amplitude_acceptance_count(k)++;
  }
};

void gp::sampleHyperParametersKthComponent(
    uword k, 
    uword n_k, 
    vec mu_tilde, 
    vec sample_mean,
    mat cov_tilde
) {
  sampleAmplitude(
    k,
    n_k,
    mu_tilde,
    sample_mean,
    cov_tilde
  );

  sampleLength(
    k,
    n_k,
    mu_tilde,
    sample_mean,
    cov_tilde
  );
  
  kernel_sub_block.slice(k) = calculateKthComponentKernelSubBlock(
    amplitude(k),
    length(k)
  );
};

double gp::noiseLogKernel(uword n_k, double noise, vec mean_vec, mat data) {
  double score = 0.0, item_score = 0.0, prior_contribution = 0.0;
  for(uword n = 0; n < n_k; n++) {
    // score += pNorm(data.row(n).t(), mean_vec, noise * I_p, true);
    item_score = 0.0;
    for(uword p = 0; p < P; p++) {
      // Normal log likelihood
      item_score -= 0.5 * std::pow(data(n, p) - mean_vec(p), 2.0);
    }
    item_score *= 1.0 / noise;
    item_score -= 0.5 * (double) P * (log(2.0 * M_PI) + log(noise));
    score += item_score;
  }
  prior_contribution = noisePriorLogDensity(noise, logNormPriorUsed); 
  score += prior_contribution;
  return score;
};

void gp::sampleNoise(uword k, uword n_k, mat component_data, double threshold) {
  bool accept = false;
  double 
      acceptance_prob = 0.0, 
      new_score = 0.0, 
      old_score = 0.0,
      new_noise = 0.0;

  new_noise = proposeNewNonNegativeValue(
    noise(k), 
    noise_proposal_window,
    use_log_norm_proposal
  );
  
  if(new_noise < threshold) {
    return;
  }
  
  new_score = noiseLogKernel(n_k, new_noise, mu.col(k), component_data);
  old_score = noiseLogKernel(n_k, noise(k), mu.col(k), component_data);
  
  acceptance_prob =  std::min(1.0, std::exp(new_score - old_score));
  accept = metropolisAcceptanceStep(acceptance_prob);
  if(accept) {
    noise(k) = new_noise;
    noise_acceptance_count(k)++;
  }
};

// === Log-likelihoods =========================================================
// The log likelihood of a item belonging to each cluster.
arma::vec gp::itemLogLikelihood(arma::vec item) {
  double exponent = 0.0;
  arma::vec ll(K), dist_to_mean(P);
  mat noise_matrix(P, P), inverse_noise_matrix(P, P);
  ll.zeros();
  dist_to_mean.zeros();
  for(uword k = 0; k < K; k++) {  
    ll(k) = logLikelihood(item, k);
  }
  return(ll);
};

// The log likelihood of a item belonging to a specific cluster.
double gp::logLikelihood(arma::vec item, arma::uword k) {
  double ll = 0.0, dist_to_mean = 0.0, exponent = 0.0, ll1, ll2;
  
  mat std_dev(P, P);
  std_dev.zeros();
  std_dev.diag().fill(noise(k));
  
  // The exponent part of the gaussian pdf
  for(uword p = 0; p < P; p++) {
    // Normal log likelihood
    ll -= 0.5 * std::pow(item(p) - mu(p, k), 2.0);
  }
  ll *= 1.0 / noise(k);
  ll -= 0.5 * (double) P * (log(2.0 * M_PI) + log(noise(k)));

  // ll1 = pNorm(item, mu.col(k), std_dev, false);
  // ll2 = gaussianLogLikelihood(item, mu.col(k), std_dev.diag());
  // 
  // Rcpp::Rcout << "\n\npNorm: " << ll1 
  //   << "\ngaussianLogLihood:" << ll2 
  //   << "\nGP implementation " << ll;
  
  return(ll);
};
