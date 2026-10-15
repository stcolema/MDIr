// jointAllocation.cpp
// =============================================================================
// Test hook for the exact joint draw of the labels of one item across views.
# include "mdi.h"

using namespace Rcpp ;
using namespace arma ;

//' @title Draw the labels of one item jointly across views
//' @description Repeated draws from the exact joint conditional of the labels of one item
//' in a block of views (see `mdiSampleJointBlock`), for checking the draw against
//' enumeration.
//' @param log_g Log of the likelihood of the item in each component of each view (K_max x
//' L); -Inf marks a component the item cannot take.
//' @param w The weights (K_max x L).
//' @param K Number of components of each view.
//' @param phi Symmetric matrix of the phis.
//' @param block The views drawn together, from 0.
//' @param current The labels of the item in every view, from 0 (those outside the block
//' are held fixed).
//' @param n_draws Number of draws.
//' @return A matrix with a row for each draw and a column for each view of the block.
// [[Rcpp::export]]
arma::umat jointAllocationDrawCpp(arma::mat log_g, arma::mat w, arma::uvec K, arma::mat phi,
                                  arma::uvec block, arma::uvec current, arma::uword n_draws) {
  arma::umat out(n_draws, block.n_elem);
  for(uword i = 0; i < n_draws; i++) {
    out.row(i) = mdiSampleJointBlock(log_g, w, K, phi, block, current, nullptr).t();
  }
  return out;
}
