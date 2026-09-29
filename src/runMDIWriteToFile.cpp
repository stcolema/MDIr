// runMDIWriteToFile.cpp
// =============================================================================
# include "runMDIWriteToFile.h"

using namespace Rcpp ;
using namespace arma ;

// Each saved sample is one binary file holding a vector with layout
//   [ labels (N * L, view by view) | weights (sum(K), view by view) | 
//     mass (L) | phis (L (L - 1) / 2) | complete likelihood | observed likelihood ]
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
    std::string save_dir,
    arma::vec prior
) {
  
  if(thin < 1) {
    Rcpp::stop("thin must be a positive integer.");
  }
  
  mdi my_mdi(Y, mixture_types, outlier_types, K, labels, fixed, prior);
  const uword L = my_mdi.L, N = my_mdi.N, LC2 = my_mdi.LC2, K_sum = accu(K);
  const std::string gen_filename = save_dir + "/MDIMcmcSample";
  
  for(uword l = 0; l < L; l++) {
    if(mixture_types[l] == 3) {
      my_mdi.mixtures[l]->density_ptr->receiveHyperParametersProposalWindows(proposal_windows[l]);
    }
  }
  
  vec save_vec(N * L + K_sum + L + LC2 + 2);
  const uword weight_offset = N * L, mass_offset = weight_offset + K_sum, 
    phi_offset = mass_offset + L, likelihood_offset = phi_offset + LC2;
  
  auto save_sample = [&](uword s) {
    uword weight_start = weight_offset;
    for(uword l = 0; l < L; l++) {
      save_vec.subvec(l * N, (l + 1) * N - 1) = conv_to<vec>::from(my_mdi.labels.col(l));
      save_vec.subvec(weight_start, weight_start + K(l) - 1) = my_mdi.w(span(0, K(l) - 1), l);
      weight_start += K(l);
    }
    save_vec.subvec(mass_offset, mass_offset + L - 1) = my_mdi.mass;
    if(LC2 > 0 && L > 1) {
      save_vec.subvec(phi_offset, phi_offset + LC2 - 1) = my_mdi.phis;
    }
    save_vec(likelihood_offset) = my_mdi.complete_likelihood;
    save_vec(likelihood_offset + 1) = my_mdi.observed_likelihood;
    
    const std::string filename = gen_filename + std::to_string(s) + ".bin";
    if(!save_vec.save(filename)) {
      Rcpp::stop("Could not write MCMC sample to '%s'.", filename);
    }
  };
  
  my_mdi.complete_likelihood = accu(my_mdi.complete_likelihood_vec);
  my_mdi.observed_likelihood = accu(my_mdi.observed_likelihood_vec);
  uword save_ind = 0;
  save_sample(save_ind);
  
  for(uword r = 0; r < R; r++) {
    Rcpp::checkUserInterrupt();
    my_mdi.sweep(r);
    if((r + 1) % thin == 0) {
      save_ind++;
      save_sample(save_ind);
    }
  }
};
