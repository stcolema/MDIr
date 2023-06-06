// runMDIWriteToFile.cpp
// =============================================================================
// included dependencies
# include "runMDIWriteToFile.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// runMDIWriteToFile function implementation
void runMDIWriteToFile(
    arma::uword R,
    arma::uword thin,
    arma::field< arma::mat > Y,
    arma::uvec K,
    arma::uvec mixture_types,
    arma::uvec outlier_types,
    arma::umat labels,
    arma::umat fixed,
    arma::field< arma::vec > proposal_windows,
    std::string save_dir
) {
  
  // Indicator if the current iteration should be recorded
  bool save_this_iteration = false;
  
  uword L = size(Y)(0),
    n_saved = floor(R / thin) + 1,
    save_ind = 0,
    N = 0,
    LC2 = 0,
    save_vec_size = 0,
    K_sum = as_scalar(sum(K)),
    curr_col = 0;
  
  uvec label_inds,
    mass_inds, 
    weight_inds, 
    phi_inds, 
    label_start, 
    col_change(4), 
    K_cum, 
    weight_start;
  col_change.zeros();
  
  vec save_vec;
  
  std::string gen_filename = save_dir + "/MDIMcmcSample", curr_filename;
  
  
  field<mat> X(L);
  for(uword l = 0; l < L; l++) {
    X(l) = Y(l);
  }
  
  mdi my_mdi(X, mixture_types, outlier_types, K, labels, fixed);
  N = my_mdi.N;
  LC2 = my_mdi.LC2;
  K_cum = cumsum(K);
  
  // The save vector contains the item labels from each dataset (N * L), the 
  // component weights (sum(K)) and masses (L), and the phi vector (LC2)
  save_vec_size = N * L + K_sum + L + LC2;
  label_inds = linspace< uvec >(curr_col, N * L, N * L);
  curr_col += my_mdi.N * L;
  col_change(1) = curr_col;
  
  weight_inds = linspace< uvec >(curr_col, curr_col + K_sum, K_sum);
  curr_col += K_sum;
  col_change(2) = curr_col;
  
  mass_inds = linspace< uvec >(curr_col, curr_col + L, L);
  curr_col += L;
  col_change(3) = curr_col;
  
  phi_inds = linspace< uvec >(curr_col, curr_col + LC2, LC2);
  
  label_start = linspace< uvec >(0, N * (L - 1), L);
  weight_start = K_cum - K(0);
  save_vec.set_size(save_vec_size);
  save_vec.zeros();
  
  for(uword l = 0; l < L; l++) {
    if(mixture_types[l] == 3) {
      my_mdi.mixtures[l]->density_ptr->receiveHyperParametersProposalWindows(proposal_windows[l]);
    }
  }
  
  // Save the initial values for each object
  for(uword l = 0; l < L; l++) {
    save_vec.subvec(label_start(l), label_start(l) + N - 1) = 
      conv_to<vec>::from(my_mdi.labels.col(l));
    
    save_vec.subvec(col_change(1) + weight_start(l), col_change(1) + K_cum(l) - 1) = my_mdi.w.col(l);
    
  }
  save_vec.subvec(col_change(2), col_change(2) + L - 1) = my_mdi.mass;
  save_vec.subvec(col_change(3), col_change(3) + LC2 - 1) = my_mdi.phis;
  
  curr_filename = gen_filename + std::to_string(save_ind) + ".bin";
  save_vec.save(curr_filename);
  Rcpp::Rcout << "\nFile name:n" << curr_filename;
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
        save_vec.subvec(label_start(l), label_start(l) + N - 1) =  
          conv_to<vec>::from(my_mdi.labels.col(l));
        save_vec.subvec(col_change(1) + weight_start(l), col_change(1) + K_cum(l) - 1) = my_mdi.w.col(l);
        
      }
      save_vec.subvec(col_change(2), col_change(2) + L - 1) = my_mdi.mass;
      save_vec.subvec(col_change(3), col_change(3) + LC2 - 1) = my_mdi.phis;
      
      curr_filename = gen_filename + std::to_string(save_ind) + ".bin";
      save_vec.save(curr_filename);
      
    }
  }
};