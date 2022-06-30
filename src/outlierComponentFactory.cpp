// outlierComponentFactory.cpp
// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "outlierComponentFactory.h"

using namespace arma ;

// =============================================================================
// virtual outlierComponentFactory class

// empty contructor
outlierComponentFactory::outlierComponentFactory() { };
outlierComponentFactory::outlierComponentFactory(const outlierComponentFactory &L) { };

std::unique_ptr<outlierComponent> outlierComponentFactory::createOutlierComponent(
    outlierType type, arma::uvec fixed, arma::mat X
) {
  // Rcpp::Rcout << "\nMake component.";
  
  switch (type) {
  case E: {
    // Rcpp::Rcout << "\nMake empty component.";
    return std::make_unique<noOutliers>(fixed, X);
  }
  case MVT: {
    // Rcpp::Rcout << "\nMake MVT component.";
    return std::make_unique<mvt>(fixed, X);
  }
  default : {
    // Rcpp::Rcout << "\nThrow an error.";
    Rcpp::Rcerr << "invalid outlier type.\n";
    throw;
  }
  }
};

