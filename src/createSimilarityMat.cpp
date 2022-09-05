// createSimilarityMat.cpp
// =============================================================================
// included dependencies
# include "createSimilarityMat.h"

// =============================================================================
// namespace
using namespace Rcpp ;
using namespace arma ;

// =============================================================================

arma::mat createSimilarityMat(arma::umat allocations){
  
  double entry = 0.0;                     // Hold current value
  uword N = allocations.n_cols;           // Number of items clustered
  uword n_iter = allocations.n_rows;      // Number of MCMC samples taken
  mat out = ones < mat > (N, N);          // Output similarity matrix 
  
  // Compare every entry to every other entry. As symmetric and diagonal is I
  // do not need to compare points with self and only need to calculate (i, j) 
  // entry
  for (arma::uword i = 0; i < N - 1; i++){ 
    for (arma::uword j = i + 1; j < N; j++){
      entry = (double)sum(allocations.col(i) == allocations.col(j)) / n_iter ;
      out(i, j) = entry;
      out(j, i) = entry;
    }
  }
  return out;
};
