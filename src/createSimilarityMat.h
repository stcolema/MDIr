// createSimilarityMat.h
// =============================================================================
// include guard
#ifndef CREATESIMILARITYMAT_H
#define CREATESIMILARITYMAT_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
using namespace arma ;

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// createSimilarityMat function header

//' @title Create Similarity Matrix
//' @description Constructs a similarity matrix of the pairwise coclustering 
//' rate.
//' @param allocations Matrix of sampled partitions. Columns correspond to 
//' items/samples being clustered, each row is a sampled partition. 
//' @return A symmetric N x N matrix (for N rows in allocations) describing 
//' the fraction of iterations for which each pairwise combination of points are
//' assigned the same label.
//' @examples
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
arma::mat createSimilarityMat(arma::umat allocations);

#endif /* CREATESIMILARITYMAT_H */
