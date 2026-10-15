// runMDIWriteToFile.h
// =============================================================================
// include guard
#ifndef RUNMDIWRITETOFILE_H
#define RUNMDIWRITETOFILE_H

// =============================================================================
// included dependencies
# include "logLikelihoods.h"
# include "genericFunctions.h"
# include "mdi.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// runMDI function header

//' @title Call Multiple Dataset Integration and Write to File
//' @description C++ function to perform MCMC sampling for MDI, writing samples
//' to file rather then saving to memory.
//' @param R The number of iterations to run for.
//' @param thin thinning factor for samples recorded.
//' @param Y The list of data matrices to perform integrative clustering upon 
//' with items to cluster in rows.
//' @param K Vector of the number of components to model in each view. This is the upper 
//' limit on the number of clusters that can be found.
//' @param mixture_types Character vector of densities used in each view
//' @param outlier_types Character vector of outlier components used in each 
//' view ('MVT' or 'None').
//' @param labels Matrix item labels to initialise from. Rows correspond to the
//' items being clustered, columns to views.
//' @param fixed Binary matrix of the items that are fixed in their initial
//' label.
//' @param proposal_windows List/field of vectors
//' @param save_dir Directory to save MCMC samples to
//' @param prior Optional vector of MDI-level prior hyperparameters (see `runMDI`).
//' @param density_prior Options of the density-level priors (see `runMDI`).
//' @return Nothing; one binary file per saved sample is written to `save_dir`.
// [[Rcpp::export]]
void runMDIWriteToFile(
     arma::uword R,
     arma::uword thin,
     arma::field<arma::mat> Y,
     arma::uvec K,
     arma::uvec mixture_types,
     arma::uvec outlier_types,
     arma::umat labels,
     arma::umat fixed,
     arma::field< arma::vec > proposal_windows,
     std::string save_dir,
     arma::vec prior,
     arma::vec density_prior
 );

#endif /* RUNMDIWRITETOFILE_H */