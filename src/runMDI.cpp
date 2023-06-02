// runMDI.cpp
// =============================================================================
// included dependencies
# include "runMDI.h"

// =============================================================================
// namespace
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
    arma::field< arma::vec > proposal_windows
) {
  
  // Indicator if the current iteration should be recorded
  bool save_this_iteration = false;
  
  uword L = size(Y)(0),
    n_saved = floor(R / thin) + 1,
    save_ind = 0,
    N = 0;
  
  field<mat> X(L);
  
  for(uword l = 0; l < L; l++) {
    X(l) = Y(l);
  }
  
  // mdiModel my_mdi(X, mixture_types, K, labels, fixed);
  mdi my_mdi(X, mixture_types, outlier_types, K, labels, fixed);
  
  for(uword l = 0; l < L; l++) {
    if(mixture_types[l] == 3) {
      my_mdi.mixtures[l]->density_ptr->receiveHyperParametersProposalWindows(proposal_windows[l]);
    }
  }
  
  N = my_mdi.N;
  
  vec likelihood_record(n_saved), evidence(n_saved - 1);
  
  mat phis_record(n_saved, my_mdi.LC2), mass_record(n_saved, L); //,
    // likelihood_record(n_saved, L);
  
  phis_record.zeros();
  mass_record.zeros();
  
  likelihood_record.zeros();
  evidence.zeros();
  
  ucube class_record(n_saved, N, L),
    outlier_record(n_saved, N, L),
    N_k_record(my_mdi.K_max, L, n_saved);
  
  class_record.zeros();
  outlier_record.zeros();
  N_k_record.zeros();
  
  cube weight_record(n_saved, my_mdi.K_max, L);
  weight_record.zeros();
  
  // field<mat> alloc(L);
  field< vec > acceptance_count(L);
  field < mat >  hyper_record(L);
  field< cube > alloc(L);
  
  // Progress p(n_saved, display_progress);
  
  for(uword l = 0; l < L; l++) {
    // alloc(l) = zeros<mat>(N, K(l));
    alloc(l) = zeros<cube>(N, K(l), n_saved);
    hyper_record(l) = zeros< mat >(n_saved, 3 * K(l));
    acceptance_count(l) = zeros< vec >(3 * K(l));
  }

  // Save the initial values for each object
  for(uword l = 0; l < L; l++) {
    // urowvec labels_l = my_mdi.labels.col(l).t();
    class_record.slice(l).row(save_ind) = my_mdi.labels.col(l).t();
    weight_record.slice(l).row(save_ind) = my_mdi.w.col(l).t();
    
    // Save the allocation probabilities
    // alloc(l) += my_mdi.mixtures[l]->alloc;
    alloc(l).slice(save_ind) = my_mdi.mixtures[l]->alloc;
    
    // Save the record of which items are considered outliers
    outlier_record.slice(l).row(save_ind) = my_mdi.mixtures[l]->outliers.t();
    
    // Save the complete likelihood
    // likelihood_record(save_ind, l) = my_mdi.mixtures[l]->complete_likelihood;
    
    
    if(mixture_types(l) == 3) {
      hyper_record(l).row(save_ind) = my_mdi.mixtures[l]->density_ptr->hypers.t();
      // acceptance_count(l) = my_mdi.mixtures[l]->density_ptr->acceptance_count.t();
    }
    
  }
  
  likelihood_record(save_ind) = my_mdi.complete_likelihood;
  mass_record.row(save_ind) = my_mdi.mass.t();
  phis_record.row(save_ind) = my_mdi.phis.t();
  N_k_record.slice(save_ind) = my_mdi.N_k;
  
  for(uword r = 0; r < R; r++) {
    
    Rcpp::checkUserInterrupt();
    
    // Should the current MCMC iteration be saved?
    save_this_iteration = ((r + 1) % thin == 0);
    my_mdi.updateNormalisingConstantNaive();
    my_mdi.sampleStrategicLatentVariable();
    my_mdi.updateMassParameters();
    my_mdi.updateWeights();
    my_mdi.updatePhis();
    for(uword l = 0; l < L; l++) {
      my_mdi.mixtures[l]->sampleParameters();
    }
    my_mdi.updateAllocation();

    // Try and swap labels within datasets to improve the correlation between 
    // clusterings across datasets
    if((r + 1) %  10 == 0) {
      my_mdi.updateLabels();
    }
    
    if( save_this_iteration ) {
      save_ind++;
      for(uword l = 0; l < L; l++) {
        class_record.slice(l).row(save_ind) = my_mdi.labels.col(l).t();
        weight_record.slice(l).row(save_ind) = my_mdi.w.col(l).t();
        
        // Save the allocation probabilities
        alloc(l).slice(save_ind) = my_mdi.mixtures[l]->alloc;
        
        // Save the record of which items are considered outliers
        outlier_record.slice(l).row(save_ind) = my_mdi.mixtures[l]->outliers.t();
        
        // Save the complete likelihood
        // Rcpp::Rcout << "\nSave model likelihood.";
        // likelihood_record(save_ind, l) = my_mdi.mixtures[l]->complete_likelihood;
        
        if(mixture_types(l) == 3) {
          hyper_record(l).row(save_ind) = my_mdi.mixtures[l]->density_ptr->hypers.t();
          // acceptance_count(l) = conv_to< vec >::from(my_mdi.mixtures[l]->density_ptr->acceptance_count).t() / r;
        }
      }
      
      evidence(save_ind - 1) = my_mdi.Z;
      likelihood_record(save_ind) = my_mdi.complete_likelihood;
      mass_record.row(save_ind) = my_mdi.mass.t();
      phis_record.row(save_ind) = my_mdi.phis.t();
      N_k_record.slice(save_ind) = my_mdi.N_k;
    }
  }
  for(uword l = 0; l < L; l++) {
    if(mixture_types(l) == 3) {
      acceptance_count(l) = conv_to< vec >::from(my_mdi.mixtures[l]->density_ptr->acceptance_count) / ( 0.2 * (double) R );
    }
  }
  
  return(
    List::create(
      Named("allocations") = class_record,
      Named("phis") = phis_record,
      Named("weights") = weight_record,
      Named("mass") = mass_record,
      Named("outliers") = outlier_record,
      Named("allocation_probabilities") = alloc,
      Named("N_k") = N_k_record,
      Named("complete_likelihood") = likelihood_record,
      Named("evidence") = evidence,
      Named("hypers") = hyper_record,
      Named("acceptance_count") = acceptance_count
    )
  );
  
}