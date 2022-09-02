# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma ;

//' @title Create Similarity Matrix
//' @description Constructs a similarity matrix of the pairwise coclustering 
//' rate.
//' @param allocations Matrix of sampled partitions. Columns correspond to 
//' items/samples being clustered, each row is a sampled partition.//' 
//' @return A symmetric n x n matrix (for n rows in cluster record) describing 
//' the fraction of iterations for which each pairwise combination of points are
//' assigned the same label.
//' @example
//' 
//' N <- 100
//' X <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
//' Y <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
//'
//' truth <- c(rep(1, N / 2), rep(2, N / 2))
//' data_modelled <- list(X, Y)
//'
//' V <- length(data_modelled)
//'
//' # This R is much too low for real applications
//' R <- 100
//' thin <- 5
//' burn <- 10
//'
//' K_max <- 10
//' K <- rep(K_max, V)
//' types <- rep("G", V)
//'
//' mcmc_out <- callMDI(data_modelled, R, thin, types, K = K)
//' createSimilarityMat(mcmc_out$allocations[, , 1])
//' @export
// [[Rcpp::export]]
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
}
