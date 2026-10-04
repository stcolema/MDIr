// splitMerge.cpp
// =============================================================================
// Test hooks for the collapsed marginal likelihoods and the split-merge move.
# include "mdi.h"

using namespace Rcpp ;
using namespace arma ;

//' @title Collapsed log marginal likelihood of a set of items (test hook)
//' @description The log of the integral of the likelihood of the items to the power
//' `beta` against the density's conjugate prior (hyperparameters at their
//' data-driven values).
//' @param X Data matrix.
//' @param K Number of components (enters the data-driven scale).
//' @param mixture_type Integer density code (0 = G, 1 = MVN, 2 = C).
//' @param density_prior Density-level prior options.
//' @param rows Zero-based indices of the items.
//' @param beta Inverse temperature.
//' @keywords internal
// [[Rcpp::export]]
double collapsedLogMarginalCpp(arma::mat X, arma::uword K, arma::uword mixture_type,
                               arma::vec density_prior, arma::uvec rows, double beta) {
  arma::uvec labels(X.n_rows, arma::fill::zeros);
  densityFactory::densityType val = static_cast<densityFactory::densityType>(mixture_type);
  std::unique_ptr<density> d = densityFactory::createDensity(val, K, labels, X, density_prior);
  if(!d->hasCollapsedMarginal()) {
    Rcpp::stop("No collapsed marginal for this density.");
  }
  collapsedStats st = d->emptyStats();
  for(uword i = 0; i < rows.n_elem; i++) {
    d->addItemToStats(st, rows(i));
  }
  return d->logMarginalLikelihood(st, beta);
}

//' @title Run the split-merge move on its own (test hook)
//' @description One view, weights held fixed at `w`, nothing else updated: the
//' move followed by the redraw of the parameters of the two components. The labels then
//' follow the collapsed target proportional to the product over items of `w[c_n]` and
//' over components of the collapsed marginal likelihood.
//' @param X Data matrix.
//' @param K Number of components.
//' @param mixture_type Integer density code (0 = G, 1 = MVN, 2 = C).
//' @param density_prior Density-level prior options.
//' @param labels Initial labels (zero-based).
//' @param fixed Indicator (0/1) of items with an observed label.
//' @param w Component weights (length K).
//' @param n_iter Number of move attempts.
//' @param beta Inverse temperature.
//' @return The labels after every attempt, and the acceptance rate.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List splitMergeOnlyCpp(arma::mat X, arma::uword K, arma::uword mixture_type,
                             arma::vec density_prior, arma::uvec labels, arma::uvec fixed,
                             arma::vec w, arma::uword n_iter, double beta) {
  const uword N = X.n_rows;
  arma::field<arma::mat> Y(1);
  Y(0) = X;
  arma::uvec K_vec(1), types(1), outliers(1);
  K_vec(0) = K;
  types(0) = mixture_type;
  outliers(0) = 0;
  arma::umat lab(N, 1), fix(N, 1);
  lab.col(0) = labels;
  fix.col(0) = fixed;
  mdi model(Y, types, outliers, K_vec, lab, fix, arma::vec(), density_prior);
  if(beta != 1.0) {
    model.setBeta(beta);
  }
  model.w.zeros();
  for(uword k = 0; k < K; k++) {
    model.w(k, 0) = w(k);
  }
  model.labels.col(0) = labels;
  model.mixtures[0]->labels = labels;
  model.mixtures[0]->density_ptr->labels = labels;
  model.refreshMembersViewL(0);
  model.setSplitMerge(1);
  arma::umat trace(n_iter, N);
  for(uword it = 0; it < n_iter; it++) {
    Rcpp::checkUserInterrupt();
    model.updateSplitMerge();
    trace.row(it) = model.labels.col(0).t();
  }
  return Rcpp::List::create(
    Named("labels") = trace,
    Named("acceptance") = (double) model.split_merge_accepts / std::max(1.0, (double) model.split_merge_attempts)
  );
}
