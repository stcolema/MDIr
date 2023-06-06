// readMCMCsamples.cpp
// =============================================================================
// included dependencies
# include "readMCMCsamples.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// readMCMCsamples function implementation
arma::mat readMCMCsamples(
    arma::uword n_samples,
    arma::uword n_params,
    std::string load_dir
) {
  std::string gen_filename = load_dir + "/MDIMcmcSample", curr_filename;
  vec read_vec(n_params);
  mat samples(n_samples, n_params);
  
  for(uword ii = 0; ii < n_samples; ii++) {
    curr_filename = gen_filename + std::to_string(ii) + ".bin";
    read_vec.load(curr_filename);
    samples.row(ii) = read_vec.t();
  }
  return samples;
};