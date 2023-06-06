// readMCMCsamples.h
// =============================================================================
// include guard
#ifndef READMCMCSAMPLES_H
#define READMCMCSAMPLES_H

// =============================================================================
# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// readMCMCsamples function header

//' @title Read MCMC Samples
//' @description C++ function to read the files saved from 
//' `runMDIWriteToFile.cpp` to disk and compile them into a matrix.
//' @param n_samples The number of MCMC samples to read in.
//' @param n_params The number of parameters recorded in each MCMC sample.
//' @param load_dir The directory containing the output of `runMDIWriteToFile`.
//' @return Matrix composed of the MCMC samples from `runMDIWriteToFile`.
// [[Rcpp::export]]
arma::mat readMCMCsamples(
    arma::uword n_samples,
    arma::uword n_params,
    std::string load_dir
);

#endif /* READMCMCSAMPLES_H */