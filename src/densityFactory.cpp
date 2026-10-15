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
  arma::vec density_prior
) {
  switch (type) {
    case G: return std::unique_ptr<gaussian>(new gaussian(K, labels, X, density_prior));
    case MVN: return std::unique_ptr<mvn>(new mvn(K, labels, X, density_prior));
    case C: return std::unique_ptr<categorical>(new categorical(K, labels, X));
    case GP: return std::unique_ptr<gp>(new gp(K, labels, X, density_prior));
  default : {
      Rcpp::stop("invalid density type.");
    }
  }
};
