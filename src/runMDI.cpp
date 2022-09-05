// runMDI.cpp
// =============================================================================
// included dependencies
# include "runMDI.h"

// =============================================================================
// runMDI function implementation
Rcpp::List runMDI(
    arma::uword R,
    arma::uword thin,
    arma::field<arma::mat> Y,
    arma::uvec K,
    arma::uvec mixture_types,
    arma::uvec outlier_types,
    arma::umat labels,
    arma::umat fixed
) {
  
  // Indicator if the current iteration should be recorded
  bool save_this_iteration = false;
  
  uword L = size(Y)(0),
    n_saved = floor(R / thin) + 1,
    save_ind = 0,
    N = 0;
  
  field<mat> X(L);
  
  // Rcpp::Rcout << "\nNumber of saved samples: " << n_saved;
  
  for(uword l = 0; l < L; l++) {
    X(l) = Y(l);
  }
  
  // mdiModel my_mdi(X, mixture_types, K, labels, fixed);
  mdi my_mdi(X, mixture_types, outlier_types, K, labels, fixed);
  
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
  field<cube> alloc(L);
  
  // Progress p(n_saved, display_progress);
  
  for(uword l = 0; l < L; l++) {
    // alloc(l) = zeros<mat>(N, K(l));
    alloc(l) = zeros<cube>(N, K(l), n_saved);
    
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
  }
  
  likelihood_record(save_ind) = my_mdi.complete_likelihood;
  
  mass_record.row(save_ind) = my_mdi.mass.t();
  
  // Rcpp::Rcout << "\n\nSave phis.";
  phis_record.row(save_ind) = my_mdi.phis.t();
  
  N_k_record.slice(save_ind) = my_mdi.N_k;
  
  for(uword r = 0; r < R; r++) {
    
    Rcpp::checkUserInterrupt();
    
    // Should the current MCMC iteration be saved?
    save_this_iteration = ((r + 1) % thin == 0);
    
    // Rcpp::Rcout << "\nUpdate normalising constant.";
    
    // Rcpp::Rcout << "\n\nNormalising constant.";
    my_mdi.updateNormalisingConstantNaive();
    
    // Rcpp::Rcout << "\nSample strategic latent variable.";
    
    // Rcpp::Rcout << "\nStrategic latent variable.";
    my_mdi.sampleStrategicLatentVariable();
    
    // Rcpp::Rcout << "\nSample component weights.";
    
    my_mdi.updateMassParameters();
    
    // Rcpp::Rcout << "\nWeights update.";
    my_mdi.updateWeights();
    
    // Rcpp::Rcout << "\nSample phis.";
    
    // Rcpp::Rcout << "\nPhis update.";
    my_mdi.updatePhis();
    
    // Rcpp::Rcout << "\nSample mixture parameters.";
    
    // Rcpp::Rcout << "\nSample mixture parameters.";
    for(uword l = 0; l < L; l++) {
      my_mdi.mixtures[l]->sampleParameters();
    }
    
    // Rcpp::Rcout << "\nUpdate allocation.";
    
    my_mdi.updateAllocation();

    // Rcpp::Rcout << "\nSwap labels.";
    
    // Try and swap labels within datasets to improve the correlation between 
    // clusterings across datasets
    if((r + 1) %  10 == 0) {
      my_mdi.updateLabels();
    }
    
    if( save_this_iteration ) {
      save_ind++;

      // Rcpp::Rcout << "\nSave objects.";
      for(uword l = 0; l < L; l++) {
        
        // urowvec labels_l = my_mdi.labels.col(l).t();
        class_record.slice(l).row(save_ind) = my_mdi.labels.col(l).t();
        
        // Rcpp::Rcout << "\nSave weights.";
        
        weight_record.slice(l).row(save_ind) = my_mdi.w.col(l).t();
        
        // Save the allocation probabilities
        // alloc(l) += my_mdi.mixtures[l]->alloc;
        // Rcpp::Rcout << "\nSave allocation probabilities.";
        
        alloc(l).slice(save_ind) = my_mdi.mixtures[l]->alloc;
        
        // Save the record of which items are considered outliers
        // Rcpp::Rcout << "\nSave outliers.";
        outlier_record.slice(l).row(save_ind) = my_mdi.mixtures[l]->outliers.t();
        
        // Save the complete likelihood
        // Rcpp::Rcout << "\nSave model likelihood.";
        // likelihood_record(save_ind, l) = my_mdi.mixtures[l]->complete_likelihood;
      }
      
      evidence(save_ind - 1) = my_mdi.Z;
      
      likelihood_record(save_ind) = my_mdi.complete_likelihood;
      
      mass_record.row(save_ind) = my_mdi.mass.t();
      
      // Rcpp::Rcout << "\n\nSave phis.";
      phis_record.row(save_ind) = my_mdi.phis.t();
      
      N_k_record.slice(save_ind) = my_mdi.N_k;
    }
    
    // p.increment(); 
    
    // Rcpp::Rcout << r << "th iteration done.\n";
    // throw;
  }
  // Rcpp::Rcout << "\nNumber of times accepted: " << my_mdi.acceptance_count << "\nPossible acceptance: " <<
  //   arma::accu(R * my_mdi.K);
  
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
      Named("evidence") = evidence
    )
  );
  
}