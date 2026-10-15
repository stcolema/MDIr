// densityFactory.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "densityFactory.h"

using namespace arma ;

// =============================================================================
// virtual densityFactory class

// empty contructor
densityFactory::densityFactory() { };
densityFactory::densityFactory(const densityFactory &L) { };

std::unique_ptr<density> densityFactory::createDensity(
  densityType type,
  arma::uword K,
  arma::uvec labels,
  arma::mat X,
  arma::vec density_prior,
  arma::vec time_points
) {
  if(!time_points.is_empty() && type != GP) {
    Rcpp::stop("`time_points` applies to Gaussian process views only.");
  }
  switch (type) {
    case G: return std::unique_ptr<gaussian>(new gaussian(K, labels, X, density_prior));
    case MVN: return std::unique_ptr<mvn>(new mvn(K, labels, X, density_prior));
    case C: return std::unique_ptr<categorical>(new categorical(K, labels, X));
    case GP: return std::unique_ptr<gp>(new gp(K, labels, X, density_prior, time_points));
  default : {
      Rcpp::stop("invalid density type.");
    }
  }
};
