// densityFactory.h
// =============================================================================
// include guard
#ifndef DENSITYFACTORY_H
#define DENSITYFACTORY_H

// =============================================================================
// included dependencies
# include <RcppArmadillo.h>
# include "density.h"
# include "gaussian.h"
# include "mvn.h"
# include "categorical.h"
# include "gp.h"

using namespace arma ;

// =============================================================================
// virtual densityFactory class

class densityFactory {
  
public:
  
  // empty contructor
  densityFactory();
  
  // destructor
  virtual ~densityFactory(){ };
  
  // copy constructor
  densityFactory(const densityFactory &L);
  
  enum densityType {
    G = 0,
    MVN = 1,
    C = 2,
    GP = 3
  };
  
  static std::unique_ptr<density> createDensity(
    densityType type,
    arma::uword K,
    arma::uvec labels,
    arma::mat X
  );
};

#endif /* DENSITYFACTORY_H */