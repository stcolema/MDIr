// runMDI.cpp
// =============================================================================
# include "runMDI.h"

using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// runMDI function implementation
Rcpp::List runMDI(
    arma::uword R,
    arma::uword thin,
    arma::field< arma::mat > Y,
    arma::uvec K,
    arma::uvec mixture_types,
    arma::uvec outlier_types,
    arma::umat labels,
    arma::umat fixed,
    arma::field< arma::vec > proposal_windows,
    bool save_parameters,
    bool save_imputed,
    arma::vec prior,
    arma::vec density_prior,
    arma::uvec save_allocation_probabilities,
    bool save_pointwise,
    bool phi_slice
) {
  
  if(thin < 1) {
    Rcpp::stop("thin must be a positive integer.");
  }
  
  const uword L = Y.n_elem, n_saved = R / thin + 1;
  uword save_ind = 0;
  
  if(save_allocation_probabilities.n_elem == 1) {
    save_allocation_probabilities = arma::uvec(L, arma::fill::value(save_allocation_probabilities(0)));
  }
  if(save_allocation_probabilities.n_elem != L) {
    Rcpp::stop("save_allocation_probabilities must have one entry per view.");
  }
  
  mdi my_mdi(Y, mixture_types, outlier_types, K, labels, fixed, prior, density_prior);
  my_mdi.phi_slice = phi_slice;
  
  for(uword l = 0; l < L; l++) {
    // Only Gaussian process views use proposal windows
    if(mixture_types[l] == 3) {
      my_mdi.mixtures[l]->density_ptr->receiveHyperParametersProposalWindows(proposal_windows[l]);
    }
  }
  
  const uword N = my_mdi.N;
  
  vec complete_likelihood_record(n_saved, arma::fill::zeros), 
    observed_likelihood_record(n_saved, arma::fill::zeros),
    joint_likelihood_record(n_saved, arma::fill::zeros),
    evidence(n_saved, arma::fill::zeros);
  mat pointwise_record;
  if(save_pointwise) {
    pointwise_record.zeros(n_saved, N);
  }
  
  mat phis_record(n_saved, my_mdi.LC2, arma::fill::zeros), 
    mass_record(n_saved, L, arma::fill::zeros),
    outlier_weight_record(n_saved, L, arma::fill::zeros);
  
  ucube class_record(n_saved, N, L), outlier_record(n_saved, N, L);
  class_record.zeros();
  outlier_record.zeros();
  
  cube weight_record(n_saved, my_mdi.K_max, L, arma::fill::zeros);
  ucube N_k_record(my_mdi.K_max, L, n_saved, arma::fill::zeros);
  
  field< cube > alloc(L);
  field< mat > hyper_record(L), parameter_record(L), imputed_record(L), pooled_record(L);
  field< umat > missing_cells(L);
  field< vec > acceptance_count(L);
  
  for(uword l = 0; l < L; l++) {
    // Only needed for semi-supervised views; N x K x draws is large
    alloc(l) = save_allocation_probabilities(l) ? zeros<cube>(N, K(l), n_saved) : zeros<cube>(0, 0, 0);
    hyper_record(l) = zeros< mat >(n_saved, 3 * K(l));
    pooled_record(l) = zeros< mat >(n_saved, my_mdi.mixtures[l]->density_ptr->pooledHyperparameters().n_elem);
    acceptance_count(l) = zeros< vec >(3 * K(l));
    
    if(save_parameters) {
      parameter_record(l) = zeros< mat >(n_saved, my_mdi.mixtures[l]->density_ptr->parameters().n_elem);
    }
    
    // The positions of missing entries (0-based row, column); fixed for the run
    const umat& has_missing = my_mdi.mixtures[l]->density_ptr->has_missing;
    uvec missing_linear = find(has_missing == 1);
    missing_cells(l) = zeros<umat>(missing_linear.n_elem, 2);
    for(uword i = 0; i < missing_linear.n_elem; i++) {
      missing_cells(l)(i, 0) = missing_linear(i) % N;
      missing_cells(l)(i, 1) = missing_linear(i) / N;
    }
    if(save_imputed) {
      imputed_record(l) = zeros< mat >(n_saved, missing_linear.n_elem);
    }
  }
  
  auto record = [&](uword s) {
    for(uword l = 0; l < L; l++) {
      class_record.slice(l).row(s) = my_mdi.labels.col(l).t();
      weight_record.slice(l).row(s) = my_mdi.w.col(l).t();
      if(save_allocation_probabilities(l)) {
        alloc(l).slice(s) = my_mdi.mixtures[l]->alloc;
      }
      outlier_record.slice(l).row(s) = my_mdi.mixtures[l]->outliers.t();
      outlier_weight_record(s, l) = my_mdi.mixtures[l]->outlierComponent_ptr->outlier_weight;
      
      const auto& density_ptr = my_mdi.mixtures[l]->density_ptr;
      if(mixture_types(l) == 3) {
        hyper_record(l).row(s) = density_ptr->hypers.t();
      }
      if(pooled_record(l).n_cols > 0) {
        pooled_record(l).row(s) = density_ptr->pooledHyperparameters().t();
      }
      if(save_parameters) {
        parameter_record(l).row(s) = density_ptr->parameters().t();
      }
      if(save_imputed) {
        for(uword i = 0; i < missing_cells(l).n_rows; i++) {
          imputed_record(l)(s, i) = density_ptr->X(missing_cells(l)(i, 0), missing_cells(l)(i, 1));
        }
      }
    }
    const arma::vec pointwise = my_mdi.pointwiseLogLikelihood();
    joint_likelihood_record(s) = accu(pointwise);
    if(save_pointwise) {
      pointwise_record.row(s) = pointwise.t();
    }
    complete_likelihood_record(s) = my_mdi.complete_likelihood;
    observed_likelihood_record(s) = my_mdi.observed_likelihood;
    evidence(s) = my_mdi.Z_start;
    mass_record.row(s) = my_mdi.mass.t();
    phis_record.row(s) = my_mdi.phis.t();
    N_k_record.slice(s) = my_mdi.N_k;
  };
  
  // The initial state
  my_mdi.complete_likelihood = accu(my_mdi.complete_likelihood_vec);
  my_mdi.observed_likelihood = accu(my_mdi.observed_likelihood_vec);
  record(save_ind);
  
  for(uword r = 0; r < R; r++) {
    Rcpp::checkUserInterrupt();
    my_mdi.sweep(r);
    if((r + 1) % thin == 0) {
      save_ind++;
      record(save_ind);
    }
  }
  
  for(uword l = 0; l < L; l++) {
    if(mixture_types(l) == 3) {
      acceptance_count(l) = conv_to< vec >::from(my_mdi.mixtures[l]->density_ptr->acceptance_count) 
        / std::max(1.0, std::ceil((double) R / 5.0));
    }
  }
  
  return(
    List::create(
      Named("allocations") = class_record,
      Named("phis") = phis_record,
      Named("weights") = weight_record,
      Named("mass") = mass_record,
      Named("outliers") = outlier_record,
      Named("outlier_weights") = outlier_weight_record,
      Named("allocation_probabilities") = alloc,
      Named("N_k") = N_k_record,
      Named("complete_likelihood") = complete_likelihood_record,
      Named("observed_likelihood") = observed_likelihood_record,
      Named("joint_likelihood") = joint_likelihood_record,
      Named("pointwise_likelihood") = pointwise_record,
      Named("evidence") = evidence,
      Named("hypers") = hyper_record,
      Named("acceptance_count") = acceptance_count,
      Named("mass_acceptance_rate") = my_mdi.mass_acceptance_count / (double) std::max<uword>(R, 1),
      Named("parameters") = parameter_record,
      Named("pooled_hyperparameters") = pooled_record,
      Named("imputed") = imputed_record,
      Named("missing_cells") = missing_cells
    )
  );
}
