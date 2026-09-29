// noOutliers.h
// =============================================================================
// include guard
#ifndef NOOUTLIERS_H
#define NOOUTLIERS_H

# include "outlierComponent.h"

using namespace arma ;

// =============================================================================
// noOutliers class of outlier component: the default, which absorbs nothing.
class noOutliers : virtual public outlierComponent {
  
public:
  
  noOutliers(
    arma::uvec _fixed, 
    arma::mat _X,
    const arma::field<arma::uvec>* miss_idx,
    const arma::field<arma::uvec>* obs_idx
  );
  
  virtual ~noOutliers() { };
  
  bool active() const override { return false; }
  double calculateItemLogLikelihood(arma::uword n) override;
  arma::vec sampleMissingValues(arma::uword n) const override;
  arma::vec simulate() const override;
  void updateWeights(const uvec& non_outliers, const uvec& outliers) override;
};

#endif /* NOOUTLIERS_H */
