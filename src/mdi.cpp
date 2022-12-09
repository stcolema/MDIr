// mdi.cpp
// =============================================================================
// included dependencies
# include <RcppParallel.h>
# include <RcppArmadillo.h>
# include <execution>
# include "mdi.h"

using namespace arma ;

// =============================================================================
// MDI class

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;


mdi::mdi(
  arma::field<arma::mat> _X,
  uvec _mixture_types,
  uvec _outlier_types,
  arma::uvec _K,
  arma::umat _labels,
  arma::umat _fixed
) {
  
  // Mixture types
  mixture_types = _mixture_types;
  
  // The type of outlier component used
  outlier_types = _outlier_types;
  
  // First allocate the inputs to their saved, class
  
  // The number of datasets modelled
  L = size(_X)(0);
  
  // The number of pairwise combinations
  if(L > 1) {
    LC2 = L * (L - 1) / 2;
  }
  
  // Used to check all datasets have matching number of rows
  uvec N_check(L);
  
  // The number of components modelled in each dataset
  K = _K;
  K_max = max(K);
  // K_prod = prod(K);
  
  // Various objects that are used in the book-keeping for MDI. The weights,
  // normalising constant and correlation coefficients all involve some pretty
  // intense book-keeping.
  K_to_the_L = std::pow(K_max, L);
  
  // The count of members in each cluster
  N_k.set_size(K_max, L);
  N_k.zeros();
  
  // // We want to track which combinations should be unweighed and by what phi.
  // // This object will be used in calculating the normalising constant (Z), the
  // // cluster weights (gammas) and the correlation coefficient (phis)
  // // along with the phi_indicator matrix.
  // comb_inds.set_size(K_to_the_L, L);
  // comb_inds.zeros();
  
  // Concentration mass for each view
  mass.set_size(L);
  mass.zeros();
  
  // A weight vector for each dataset. Note that for ease of manipulations we
  // are using K_max and not K(l); this is to avoid ragged fields of vectors.
  w.set_size(K_max, L);
  w.zeros();
  
  // The ``correlation'' coefficient between each dataset clustering. Use the
  // size that best corresponds to phi_indicator. Possibly a row vector would
  // be easier.
  phis.set_size(LC2);
  phis.zeros();
  
  // Initialise the different objects used to interact with the phis
  initialisePhis();
  
  // The initial labels
  labels = _labels;
  
  // The data
  X = _X;
  
  // Iterate over the number of datasets (currently must be 2)
  for(uword l = 0; l < L; l++) {
    N_check(l) = _X(l).n_rows;
  }
  
  if(any(N_check != N_check(0))) {
    Rcpp::Rcerr << "\n\nDatasets not matching in number of rows.\n";
    throw;
  }
  
  // The number of samples in each dataset
  N = N_check(0);
  
  N_inds = linspace< uvec >(0, N - 1, N);
  L_inds = linspace< uvec >(0, L - 1, L);
  
  // Quantities used in calculating phi shape
  N_ones = regspace(0, N);
  N_log_factorial_vec = log(N_ones);
  N_log_factorial_vec(0) = 0.0;
  N_log_factorial_vec = cumsum(N_log_factorial_vec);

  // The members of each cluster across datasets. Each slice is a binary matrix
  // of the members of the kth class across the datasets.
  members.set_size(N, K_max, L);
  members.zeros();
  
  
  // These are used in t-adjusted mixtures. In all other models they should
  // never be changed.
  non_outliers.set_size(N, L);
  non_outliers.ones();
  
  // The observed labels
  fixed = _fixed;
  
  // The number of fixed components (i.e. non-symbolic labels) and unfixed
  // components. Used when aligning components across views.
  K_fixed.set_size(L);
  K_unfixed.set_size(L);
  
  uvec fixed_l, fixed_labels, labels_l, fixed_components;
  
  for(uword l = 0; l < L; l++){
    
    fixed_l = find(fixed.col(l) == 1);
    labels_l = labels.col(l);
    
    
    // fixed_ind(l) = fixed_l;
    fixed_labels = labels_l(fixed_l);
    fixed_components = unique(fixed_labels);
    K_fixed(l) = fixed_components.n_elem;
    K_unfixed(l) = K(l) - K_fixed(l);
    
  }
  
  // Rcpp::Rcout << "K_fixed:\n" << K_fixed.t();
  // Rcpp::Rcout << "K_fixed:\n" << K_unfixed.t();
  
  
  complete_likelihood_vec = zeros< vec >(L);
  
  // Declare the view-specific mixtures
  initialiseMDI();
};


void mdi::initialisePhis() {
  
  uword col_ind = 0;
  
  // // We want to track which combinations should be unweighed and by what phi.
  // // This object will be used in calculating the normalising constant (Z), the
  // // cluster weights (gammas) and the correlation coefficient (phis)
  // // along with the phi_indicator matrix.
  // comb_inds.set_size(K_to_the_L, L);
  // comb_inds.zeros();
  
  // The gammas combine in a form like (letting g denote gamma)
  // g_{11} g_{12} ... g_{1L}
  //          .
  //          .
  //          .
  // g_{K1} g_{12} ... g_{1L}
  // g_{11} g_{22} ... g_{1L}
  //          .
  //          .
  //          .
  // g_{K1} g_{K2} ... g_{KL}
  // 
  // Hence our matrix is initially of size K^L \times L
  one_to_K = linspace<uvec>(0, K_max - 1, K_max);
  one_to_L = linspace<uvec>(0, L - 1, L);
  
  // // The matrix used to construct the rate for sampling the cluster weights
  // KL_powers.set_size(L);
  // for(uword l = 0; l < L; l++) {
  //   KL_powers(l) = std::pow(K_max, l);
  //   for(uword i = 0; i < K_to_the_L; i++){
  //     
  //     // We want the various combinations of the different gammas / cluster weights.
  //     // This format makes it relatively easy to figure out the upscaling too
  //     // (see phi_indicator below).
  //     k = ((double) i / (double) KL_powers(l));
  //     k = k % K_max;
  //     comb_inds(i, l) = k;
  //   }
  // }
  // 
  // // Drop any rows that contain weights for clusters that shouldn't be 
  // // modelled (i.e. if unique K_l are used)
  // for(uword l = 0; l < L; l++) {
  //   rows_to_shed = find(comb_inds.col(l) >= K(l));
  //   comb_inds.shed_rows(rows_to_shed);
  // }
  // 
  // // The final number of combinations
  // n_combinations = comb_inds.n_rows;
  // 
  // // Now construct a matrix to record which phis are upweighing which weight
  // // products, via an indicator matrix. This matrix has a column for each phi
  // // (ncol = LC2) and a row for each combination (nrow = n_combinations).
  // // This is working with the combi_inds object above. That indicated which 
  // // weights to use, this indicates the corresponding up weights (e.g.,
  // // we'd expect the first row to be all ones as all weights are for the first
  // // component in each dataset, similarly for the last row).
  // phi_indicator.set_size(n_combinations, LC2);
  // phi_indicator.zeros();
  
  // Rcpp::Rcout << "\n\nPhi indicator declared.";
  
  // Map between a dataset pair and the column index. This will be a lower
  // triangular matrix of unsigned ints
  phi_map.set_size(L, L);
  phi_map.fill(100);
  // phi_map.diag().fill(100);
  
  // Column index is, for some pair of datasets l and m,
  // \sum_{i = 0}^{l} (L - i - 1) + (m - l - 1). As this is awkward, let's
  // just use col_ind++.
  col_ind = 0;
  
  // Iterate over dataset pairings
  for(uword l = 0; l < (L - 1); l++) {
    for(uword m = l + 1; m < L; m++) {
      // Record which column index maps to which phi
      phi_map(m, l) = col_ind;
      phi_map(l, m) = col_ind;
      
      // The column index is awkward, this is the  easiest solution
      col_ind++;
      
    }
  }
  
  // // We use the transpose a surprising amount to ensure correct typing
  // phi_indicator_t = phi_indicator.t();
  // 
  // // And we often multiply by doubles, so convert to a matrix of doubles.
  // phi_indicator_t_mat = conv_to<mat>::from(phi_indicator_t);
  
};

// This wrapper to declare the mixtures at the dataset level is kept as its
// own separate function to make a semi-supervised class easier
void mdi::initialiseMixtures() {

  // Initialise the collection of mixtures
  mixtures.reserve(L);
  for(uword l = 0; l < L; l++) {
    mixtures.push_back(
      std::make_unique<mixtureModel>(
        mixture_types(l),
        outlier_types(l),
        K(l),
        labels.col(l),
        fixed.col(l),
        X(l)
      )
    );
    
    // We have to pass this back up
    non_outliers.col(l) = mixtures[l]->non_outliers;
  }
};


// double mdi::calcWeightRate(uword lstar, uword kstar) {
//   // The rate for the (k, l)th cluster weight is the sum across all of the clusters
//   // in each dataset (excluding the lth) of the product of the kth cluster
//   // weight in each of the L datasets (excluding the lth) upweigthed by
//   // the pairwise correlation across all LC2 dataset combinations.
//   
//   // The rate we return.
//   double rate = 0.0;
//   
//   // A chunk of these objects are used as our combinations matrix includes the
//   // information for k != kstar and l != lstar, so we can shed some data.
//   // Probably some of these can be dropped and calculated once at the class // [[Rcpp::export]]
//   // level, but to allow K_l != K_m I suspect it will be different for each
//   // dataset (I haven't done the maths though) and thus the objects might be
//   // awkward
//   uword n_used = 0;
//   uvec relevant_inds;
//   vec weight_products, phi_products;
//   umat relevant_combinations;
//   mat relevant_phis, phi_prod_mat;
//   
//   relevant_inds = find(comb_inds.col(lstar) == kstar);
//   relevant_combinations = comb_inds.rows(relevant_inds);
//   relevant_phis =  phi_indicator_t_mat.cols(relevant_inds);
//   n_used = relevant_combinations.n_rows;
//   
//   weight_products.ones(n_used);
//   phi_products.ones(n_used);
//   
//   // The phi products (should be a matrix of 0's and phis)
//   phi_prod_mat = relevant_phis.each_col() % phis;
//   
//   // Add 1 to each entry to have the object ready to be multiplied
//   phi_prod_mat++;
//   
//   // Vector of the products, this should have the \prod (1 + \phi_{lm} ind(k_l = k_m))
//   // ready to multiply by the weight products
//   phi_products = prod(phi_prod_mat, 0).t();
//   
//   vec w_l(K_max);
//   for(uword l = 0; l < L; l++) {
//     if(l != lstar){
//       w_l = w.col(l);
//       weight_products = weight_products % w_l.elem(relevant_combinations.col(l));
//     }
//   }
//   
//   // The rate for the gammas
//   rate = v * accu(weight_products % phi_products);
//   
//   return rate;
// };

double mdi::calcWeightRateNaiveSingleIteration(uword kstar, 
                                                       uword lstar, 
                                                       uvec current_ks) {
  double output = 1.0;
  for(uword l = 0; l < L; l++) {
    if(l != lstar) {
      output *= w(current_ks(l), l);
    }
  }
  for(uword l = 0; l < L - 1; l++) {
    for(uword m = l + 1; m < L; m++) {
      output *= (1.0 + phis(phi_map(l, m)) * (current_ks(l) == current_ks(m)));
    }
  }
  return output;
};

double mdi::calcWeightRateNaive(uword kstar, uword lstar) {
  uword K_comb = 0, m = 0;
  double rate = 0.0;
  uvec weight_ind(L), K_cumprod(L - 1), for_loop_inds(L), K_rel(L), K_rel_cum(L - 1);
  
  // We begin at the 0th component for each view and the kth for the view being 
  // updated
  weight_ind.zeros();
  weight_ind(lstar) = kstar;
  
  // The number of components in each view bar the vth
  K_rel.zeros();
  K_rel_cum.ones();
  K_rel = K;
  K_rel.shed_row(lstar);
  
  // We will use this in the for loop to control accessing and updating of 
  // weights
  for_loop_inds = regspace<uvec>(0,  1,  L - 1);
  for_loop_inds.shed_row(lstar);
  
  // We have shed an entry (leaving L-1) and we do not need the final entry for 
  // the cumulative check as it should fill once
  if(L > 2) {
    K_rel_cum(span(1, L - 2)) = K_rel(span(0, L - 3));
  }
 
  // The cumulative product of the number of components in each view is used to
  // control updating weight indices
  K_cumprod = cumprod(K_rel_cum);
  
  // This is the number of summations to perform
  K_comb = prod(K_rel);
  
  // Rcpp::Rcout << "\n\nWEIGHT RATE\nView: " << l << "\nComponent: " << k << "\n";
  for(uword ii = 0; ii < K_comb; ii++) {
    // Rcpp::Rcout << "\n\ni: " << ii;
    // Rcpp::Rcout << "\nWeight indices:\n " << weight_ind.t();
    rate += calcWeightRateNaiveSingleIteration(kstar, lstar, weight_ind);
    for(uword jj = 0; jj < (L - 1); jj++) {
      // Which view is actually being updated (skipping the vth)
      m = for_loop_inds(jj);
      
      // If in first view we always update the weight index. Otherwise we check 
      // if we are a multiple of the cumulative product of the number of 
      // components in the preceding views and update the index if so.
      if(jj == 0) {
        weight_ind(m)++;
      } else {
        if((((ii + 1) % K_cumprod(jj)) == 0) && (ii != 0)) {
          weight_ind(m)++;
        }
      }
      if(weight_ind(m) == K(m)) {
        weight_ind(m) = 0;
      }
    }
  }
  rate *= v;
  return rate;
};

double mdi::calcPhiRateNaiveSingleIteration(uword view_i, uword view_j, uvec current_ks) {
  double output = 1.0;
  for(uword l = 0; l < L; l++) {
    output *= w(current_ks(l), l);
  }
  for(uword l = 0; l < L - 1; l++) {
    for(uword m = l + 1; m < L; m++) {
      if(m != view_j) {
        output *= (1.0 + phis(phi_map(m, l)) * (current_ks(l) == current_ks(m)));
      }
    }
  }
  for(uword l = 0; l < view_j - 1; l++) {
    output *= (1.0 + phis(phi_map(l, view_j)) * (current_ks(l) == current_ks(view_j)));
  }
  return output;
};

double mdi::calcPhiRateNaive(uword view_i, uword view_j) {
  bool shed_j = K(view_i) < K(view_j), weight_updated = false;
  uword K_comb = prod(K), l = 0;
  double rate = 0.0;
  uvec weight_ind(L), K_cumprod(L - 1), backwards_inds(L), for_loop_inds(L), K_rel(L), K_rel_cum(L - 1);
  weight_ind.zeros();
  
  K_rel.zeros();
  K_rel_cum.ones();
  K_rel = K;
  
  for_loop_inds = regspace<uvec>(0,  1,  L - 1);
  
  if(shed_j) {
    K_rel.shed_row(view_j);
    for_loop_inds.shed_row(view_j);
  } else {
    K_rel.shed_row(view_i);
    for_loop_inds.shed_row(view_i);
  }
  
  // We have shed an entry and we do not need the final entry for the cumulative check
  // K_rel_cum(span(1, L - 2)) = K_rel(span(0, L - 3));
  if(L > 2) {
    K_rel_cum(span(1, L - 2)) = K_rel(span(0, L - 3));
  }
  K_cumprod = cumprod(K_rel_cum);
  K_comb = prod(K_rel);
  
  // Rcpp::Rcout << "\n\nPHI RATE\nview i: " << view_i << "\nview j: " << view_j << "\n";
  for(uword ii = 0; ii < K_comb; ii++) {
    // Rcpp::Rcout << "\nWeight indices:\n" << weight_ind.t();
    rate += calcPhiRateNaiveSingleIteration(view_i, view_j, weight_ind);
    
    // We have to hold the index for view_i and view_j the same
    for(uword jj = 0; jj < (L - 1); jj++) {
      weight_updated = false;
      l = for_loop_inds(jj);
      // Rcpp::Rcout << "\nl: " << l;
      if(jj == 0) {
        weight_ind(l)++;
        weight_updated = true;
      } else {
        // if((ii % K_cumprod(jj) == 0) && (ii != 0)) {
        if((((ii + 1) % K_cumprod(jj)) == 0) && (ii != 0)) {
          weight_ind(l)++;
          weight_updated = true;
        }
      }
      if(weight_updated && shed_j && (l == view_i)) {
        weight_ind(view_j)++;
      } 
      if(weight_updated &&  (! shed_j) && (l == view_j)) {
        weight_ind(view_i)++;
      }
      if(weight_ind(l) == K(l)) {
        weight_ind(l) = 0;
        if(shed_j && (l == view_i)) {
          weight_ind(view_j) = 0;
        } 
        if(! shed_j && (l == view_j)) {
          weight_ind(view_i) = 0;
        }
      }
    }
  }
  rate *= v;
  return rate;
};

// // The rate for the phi coefficient between the lth and mth datasets.
// double mdi::calcPhiRate(uword lstar, uword mstar) {
//   
//   // The rate we return.
//   double rate = 0.0;
//   vec w_l;
//   
//   // A chunk of these objects are used as our combinations matrix includes the
//   // information for l = lstar and l = mstar, so we can shed some data.
//   uword n_used;
//   uvec relevant_inds;
//   vec weight_products, phi_products, relevant_phis;
//   umat relevant_combinations;
//   mat relevant_phi_inidicators, phi_prod_mat;
//   
//   relevant_inds = find(comb_inds.col(lstar) == comb_inds.col(mstar));
//   relevant_combinations = comb_inds.rows(relevant_inds);
// 
//   // We only need the relevant phi indicators
//   relevant_phi_inidicators = phi_indicator_t_mat.cols(relevant_inds);
//  
//   // Drop phi_{lstar, mstar} from both the indicator matrix and the phis vector
//   relevant_phi_inidicators.shed_row(phi_map(mstar, lstar));
//   relevant_phis = phis;
//   relevant_phis.shed_row(phi_map(mstar, lstar));
//   
//   n_used = relevant_combinations.n_rows;
//   
//   weight_products.ones(n_used);
//   phi_products.ones(n_used);
//   
//   // The phi products (should be a matrix of 0's and phis)
//   phi_prod_mat = relevant_phi_inidicators.each_col() % relevant_phis;
//   
//   // Add 1 to each entry to have the object ready to be multiplied
//   phi_prod_mat++;
//   
//   // Vector of the products, this should have the \prod (1 + \phi_{lm} ind(k_l = k_m))
//   // ready to multiply by the weight products
//   phi_products = prod(phi_prod_mat, 0).t();
// 
//   for(uword l = 0; l < L; l++) {
//     if(l != lstar){
//       w_l = w.col(l);
//       weight_products = weight_products % w_l.elem(relevant_combinations.col(l));
//     }
//   }
//   
//   // The rate for the gammas
//   rate = v * accu(weight_products % phi_products);
//   
//   return rate;
// };

void mdi::updateWeightsViewL(uword l) {
  
  double shape = 0.0, rate = 0.0, posterior_shape = 0.0, posterior_rate = 0.0;
  uvec members_lk(N);
  
  // Rcpp::Rcout << "\n\nl: " << l;
  for(uword k = 0; k < K(l); k++) {
    
    // Find how many labels have the value of k. We used to consider which
    // were outliers and which were not, but the non-outliers still 
    // contribute to the component weight, but not to the component parameters
    // and we use ot hand this down to the local mixture, mistakenly using 
    // the same value for N_k for the component parameters and the weights.
    members_lk = 1 * (labels.col(l) == k);
    members.slice(l).col(k) = members_lk;
    N_k(k, l) = accu(members_lk);
    
    // The hyperparameters
    shape = 1.0 + N_k(k, l);
    // rate = calcWeightRate(l, k);
    if(L > 1) {
      rate = calcWeightRateNaive(k, l);
    } else {
      rate = 1.0;
    }

    // Sample a new weight
    posterior_shape = (mass(l) / (double) K(l)) + shape;
    posterior_rate = w_rate_prior + rate;
    w(k, l) = rGamma(posterior_shape, posterior_rate);
    
    // Pass the allocation count down to the mixture
    // (this is used in the parameter sampling)
    mixtures[l]->members.col(k) = members_lk;
  }
  mixtures[l]->N_k = N_k(span(0, K(l) - 1), l);
}

// Update the cluster weights
void mdi::updateWeights() {
  
  // for(uword l = 0; l < L; l++) {
  //   updateWeightsViewL(l);
  // }
  std::for_each(
    std::execution::par,
    L_inds.begin(),
    L_inds.end(),
    [&](uword l) {
      updateWeightsViewL(l);
    }
  );
};

double mdi::samplePhiShape(arma::uword l, arma::uword m, double rate) {
  int N_lm = 0;
  double shape = 0.0,
    u = 0.0;
  
  uvec rel_inds_l(N), rel_inds_m(N);
  vec log_weights, weights;
  
  // rel_inds_l = labels.col(l) % non_outliers.col(l);
  // rel_inds_m = labels.col(m) % non_outliers.col(m);
  // 
  // N_lm = accu(rel_inds_l == rel_inds_m);
  
  N_lm = accu(labels.col(l) == labels.col(m));
  weights = zeros<vec>(N_lm + 1);
  log_weights = calculatePhiShapeMixtureWeights(N_lm, rate);
  
  // Normalise the weights
  weights = exp(log_weights - max(log_weights));
  weights = weights / accu(weights);
  
  // Prediction and update
  u = randu<double>( );
  shape = sum(u > cumsum(weights)) ;
  return shape; 
}

void mdi::averagePhiUpdate(arma::uword l, arma::uword m, double rate) {
  int N_lm = 0;
  uvec rel_inds_l(N), rel_inds_m(N);
  vec log_weights, weights, phis_vec;
  
  N_lm = accu(labels.col(l) == labels.col(m));
  weights = zeros<vec>(N_lm + 1);
  log_weights = zeros<vec>(N_lm + 1);
  phis_vec = zeros<vec>(N_lm + 1);
  log_weights = calculatePhiShapeMixtureWeights(N_lm, rate);
  
  // Normalise the weights
  weights = exp(log_weights - max(log_weights));
  weights = weights / accu(weights);
  
  for(uword ii = 0; ii < (N_lm + 1); ii++) {
    phis_vec(ii) = weights(ii) * rGamma(ii + phi_shape_prior, rate + phi_rate_prior);
  }
  phis(phi_map(m, l)) = accu(phis_vec);
}

arma::vec mdi::calculatePhiShapeMixtureWeights(
    int N_lm, 
    double rate
) {
  
  double r_alpha_gamma_function = 0.0,
    beta_part = 0.0,
    log_n_choose_r = 0.0;

  vec N_lm_ones(N_lm + 1), 
    N_lm_vec(N_lm + 1),
    log_weights(N_lm + 1);
  log_weights.zeros();
  
  // N_lm_ones = regspace(0,  N_lm);
  // r_log_factorial_vec = N_log_factorial_vec.subvec(0, N_lm);
  for(int r = 0; r < (N_lm + 1); r++) {
    log_n_choose_r = logChoose(N_lm, r);
    r_alpha_gamma_function = lgamma(r + phi_shape_prior);
    beta_part = ((double) r + phi_shape_prior) * std::log(rate + phi_rate_prior);
    log_weights(r) = log_n_choose_r + r_alpha_gamma_function - beta_part;
  }
  return log_weights;
};

void mdi::updatePhis() {
  if(L == 1) {
    return;
  }
  double shape = 0.0, rate = 0.0, posterior_shape = 0.0, posterior_rate = 0.0;
  for(uword l = 0; l < (L - 1); l++) {
    for(uword m = l + 1; m < L; m++) {
      // Find the parameters based on the likelihood
      // rate = calcPhiRate(l, m);
      // shape = 1 + accu(labels.col(l) == labels.col(m)); // original implementationshape = 1 + accu(labels.col(l) == labels.col(m));
      rate = calcPhiRateNaive(l, m);
      shape = samplePhiShape(l, m, rate);
      posterior_shape = phi_shape_prior + shape;
      posterior_rate = phi_rate_prior + rate;
      phis(phi_map(m, l)) = rGamma(posterior_shape, posterior_rate);
      // averagePhiUpdate(l, m, rate); // smoother update
    }
  }
};

// // Update the context similarity parameters
// void updatePhis() {
//   
//   // Rcpp::Rcout << "\n\nPhis before update:\n" << phis;
//   uword N_lm = 0;
//   double shape = 0.0, rate = 0.0;
//   uvec rel_inds_l(N), rel_inds_m(N);
//   
//   for(uword l = 0; l < (L - 1); l++) {
//     for(uword m = l + 1; m < L; m++) {
//       
//       // I considered excluding the outliers from the phi calculation, but to 
//       // be consistent we would also have to exclude them from component 
//       // weights.
//       rel_inds_l = labels.col(l); // % non_outliers.col(l);
//       rel_inds_m = labels.col(m); // % non_outliers.col(m);
//       
//       N_lm = accu(rel_inds_l == rel_inds_m);
//       shape = 1 + N_lm;
//       
//       // shape = 1 + accu(labels.col(l) == labels.col(m));
//       rate = calcPhiRate(l, m);
//       
//       // Rcpp::Rcout << "\n\nShape:" << shape;
//       // Rcpp::Rcout << "\nRate:" << rate;
//       
//       phis(phi_map(m, l)) = randg(distr_param(
//         phi_shape_prior + shape,
//         1.0 / (phi_rate_prior + rate)
//       )
//       );
//     }
//   }
//   
// };

// // Updates the normalising constant for the posterior
// void mdi::updateNormalisingConst() {
//   
//   vec w_l;
//   
//   // A chunk of these objects are used as our combinations matrix includes the
//   // information for l = lstar and l = mstar, so we can shed some data.
//   vec weight_products, phi_products;
//   mat phi_prod_mat;
//   
//   weight_products.ones(n_combinations);
//   phi_products.ones(n_combinations);
//   
//   // The phi products (should be a matrix of 0's and phis)
//   phi_prod_mat = phi_indicator_t_mat.each_col() % phis;
//   
//   // Add 1 to each entry to have the object ready to be multiplied
//   phi_prod_mat++;
//   
//   // Vector of the products, this should have the \prod (1 + \phi_{lm} ind(k_l = k_m))
//   // ready to multiply by the weight products
//   phi_products = prod(phi_prod_mat, 0).t();
//   
//   
//   for(uword l = 0; l < L; l++) {
//     w_l = w.col(l);
//     
//     weight_products = weight_products % w_l.elem(comb_inds.col(l));
//   }
//   
//   // The rate for the gammas
//   Z = accu(weight_products % phi_products);
// };

double mdi::calcNormalisingConstNaiveSingleIteration(uvec current_ks) {
  double iter_value = 1.0, same_label = 0.0; // log_iter_value = 0.0, 
  
  for(uword jj = 0; jj < L; jj++) {
    iter_value *= w(current_ks(jj), jj);
    // log_iter_value += log(w(current_ks(jj), jj));
  }
  for(uword l = 0; l < L - 1; l++) {
    for(uword m = l + 1; m < L; m++) {
      same_label = 1.0 * (current_ks(l) == current_ks(m));
      iter_value *= 1.0 + phis(phi_map(m, l)) * same_label;
      // log_iter_value += log(1.0 + phis(phi_map(m, l)) * same_label);
    }
  }
  return iter_value;
};

void mdi::updateNormalisingConstantNaive() {
  uword K_comb = 0;
  uvec weight_ind(L), K_cumprod(L), K_rel(L), K_rel_cum(L);
  weight_ind.zeros();
  K_rel.zeros();
  K_rel_cum.ones();
  K_rel = K;
  
  // We do not need the final entry for the cumulative check
  if(L > 1) {
    K_rel_cum(span(1, L - 1)) = K_rel(span(0, L - 2));
  }
  K_cumprod = cumprod(K_rel_cum);
  K_comb = prod(K_rel);
  Z = 0.0;
  
  // Rcpp::Rcout << "\n\nNORMALISING CONSTANT\n";
  for(uword ii = 0; ii < K_comb; ii++) {
    // Rcpp::Rcout << "\n\ni: " << ii;
    // Rcpp::Rcout << "\nWeight indices:\n" << weight_ind.t();
    Z += calcNormalisingConstNaiveSingleIteration(weight_ind);
    for(uword l = 0; l < L; l++) {
      if(l == 0) {
        weight_ind(l)++;
      } else {
        if((((ii + 1) % K_cumprod(l)) == 0) && (ii != 0)) {
          weight_ind(l)++;
        }
      }
      if(weight_ind(l) == K(l)) {
        weight_ind(l) = 0;
      }
    }
  }
};

void mdi::sampleStrategicLatentVariable() {
  v = rGamma(N, Z);
};

void mdi::sampleFromPriors() {
  // Sample from the prior distribution for the phis and weights
  sampleFromGlobalPriors();
  
  // Sample from the prior distribution for the view-specific mixtures
  sampleFromLocalPriors();
};

void mdi::sampleFromLocalPriors() {
  for(uword l = 0; l < L; l++) {
    mixtures[l]->sampleFromPriors();
  }
};

vec mdi::samplePhiPrior(uword n_phis) {
  return rGamma(n_phis, phi_shape_prior , phi_rate_prior);
};

double mdi::sampleWeightPrior(uword l) {
  return rGamma(mass(l) / (double) K(l) , w_rate_prior);
}

vec mdi::sampleMassPrior() {
  return rGamma(L, mass_shape_prior, mass_rate_prior);
}

void mdi::sampleFromGlobalPriors() {
  mass = sampleMassPrior();
  if(L > 1) {
    phis = samplePhiPrior(LC2);
    // phis = randg(LC2, distr_param(2.0 , 1.0 / 2));
  } else {
    phis.zeros();
  }
  
  for(uword l = 0; l < L; l++) {
    for(uword k = 0; k < K(l); k++) {
      w(k, l) = sampleWeightPrior(l); // randg(distr_param(1.0 / (double)K(l), 1.0));
    }
  }
};

void mdi::updateMassParameters() {
  for(uword l = 0; l < L; l++) {
  // std::for_each(
  //   std::execution::par,
  //   L_inds.begin(),
  //   L_inds.end(),
  //   [&](uword l) {
      updateMassParameterViewL(l);
    }
  // );
}

void mdi::updateMassParameterViewL(uword lstar) {
  bool accepted = false;
  double
    current_mass = 0.0, 
    proposed_mass = 0.0, 
    
    cur_log_likelihood = 0.0,
    cur_log_prior = 0.0,
    
    new_log_likelihood = 0.0,
    new_log_prior = 0.0,
    acceptance_ratio = 0.0;
  
  vec current_weights(K(lstar));
  
  current_weights = w(span(0, K(lstar) - 1), lstar);
  current_mass = mass(lstar);
  
  // Log likelihood and log prior density
  cur_log_likelihood = gammaLogLikelihood(
    current_weights, 
    current_mass / (double) K(lstar), 
    1
  );
  
  cur_log_prior = gammaLogLikelihood(
    current_mass, 
    mass_shape_prior, 
    mass_rate_prior
  );
  
  proposed_mass = proposeNewNonNegativeValue(current_mass,
    mass_proposal_window, 
    use_log_norm_proposal
  );
  
  if(proposed_mass <= 0.0) {
    acceptance_ratio = 0.0;
  } else {
    new_log_likelihood = gammaLogLikelihood(
      current_weights, 
      proposed_mass / (double) K(lstar), 
      1.0
    );
    
    new_log_prior = gammaLogLikelihood(
      proposed_mass,  
      mass_shape_prior, 
      mass_rate_prior
    );
    
    acceptance_ratio = exp(
      new_log_likelihood
      + new_log_prior 
      - cur_log_likelihood 
      - cur_log_prior
    );
  }
  accepted = metropolisAcceptanceStep(acceptance_ratio);
  if(accepted) {
    mass(lstar) = proposed_mass;
  }
};

mat mdi::calculateUpweights(uword lstar) {
  double same_label = 0.0;
  mat log_upweights(K(lstar), N);
  log_upweights.zeros();
  // Upweights are (1 + \phi)
  for(uword m = 0; m < L; m++) {
    if(m != lstar){
      for(uword k = 0; k < K(lstar); k++) {  
        for(uword n = 0; n < N; n++ ) {
          same_label = 1.0 * (double) (labels(n, m) == k);
          log_upweights(k, n) += log(1.0 + phis(phi_map(m, lstar)) * same_label);
        }
      }
    }
  }
  return log_upweights;
};

void mdi::initialiseDatasetL(uword l) {
  vec log_weights(K(l));
  mat log_upweights(K(l), N);
  log_upweights = calculateUpweights(l);
  log_weights = log(w(span(0, K(l) - 1), l));
  mixtures[l]->initialiseMixture(log_weights, log_upweights);
  labels.col(l) = mixtures[l]->labels;
  non_outliers.col(l) = mixtures[l]->non_outliers;
}

void mdi::initialiseMDI() {
  initialiseMixtures();
  sampleFromPriors();
  for(uword l = 0; l < L; l++) {
    initialiseDatasetL(l);
  }
};

void mdi::updateAllocationViewL(uword l) {
  
  vec log_weights(K(l));
  mat log_upweights(K(l), N);
  log_weights = log(w(span(0, K(l) - 1), l));
  log_upweights = calculateUpweights(l);
  
  // Update the allocation within the mixture using MDI level weights and phis
  mixtures[l]->updateAllocation(log_weights, log_upweights);
  
  // Pass the new labels from the mixture level back to the MDI level.
  labels.col(l) = mixtures[l]->labels;
  non_outliers.col(l) = mixtures[l]->non_outliers;
  complete_likelihood_vec(l) = mixtures[l]->complete_likelihood;
}

void mdi::updateAllocation() {

  complete_likelihood = 0.0;

  std::for_each(
    std::execution::par,
    L_inds.begin(),
    L_inds.end(),
    [&](uword l) {
      updateAllocationViewL(l);
    }
  );
  complete_likelihood = accu(complete_likelihood_vec);
};


// This is used to consider possible label swaps
double mdi::sampleLabel(arma::uword kstar, arma::vec K_inv_cum) {
  
  // Need to account for the fixed labels
  // Need the non-fixed classes (e.g., we now need, K, K_fixed and K_unfixed)
  
  // Select another label randomly
  double u = randu();
  uword k_prime = sum(u > K_inv_cum);
  
  // If it is >= than the current label under consideration, add one
  if(k_prime >= kstar) {
    k_prime++;
  }
  return k_prime;
}

double mdi::calcScore(arma::uword lstar, arma::umat c) {
  bool not_current_context = true;
  double score = 0.0; // , same_label = 0.0;
  uvec agreeing_labels(N);
  vec same_label(N);
  agreeing_labels.zeros();
  same_label.zeros();
  
  for(uword m = 0; m < L; m++) {
    // Skip the current context (the lth context)
    not_current_context = (m != lstar);
    if(not_current_context) {
      same_label = conv_to<vec>::from(c.col(m) == c.col(lstar));
      score += accu(log(1.0 + phis(phi_map(m, lstar)) * same_label));
      // for(uword n = 0; n < N; n++) {
      //   same_label = 1.0 * (c(n, m) == c(n, lstar));
      //   upweight = phis(phi_map(m, lstar)) * same_label;
      //   score += log(1.0 + upweight);
      // }
    }
  }
  return score;
}

// Check if labels should be swapped to improve correlation of clustering
// across datasets via random sampling.
arma::umat mdi::swapLabels(arma::uword lstar, arma::uword kstar, arma::uword k_prime) {
  
  // The labels in the current context, which will be changed
  uvec loc_labs = labels.col(lstar), 
    
    // The indices for the clusters labelled k and k prime
    cluster_k,
    cluster_k_prime;
  
  // The labels for the other contexts
  umat dummy_labels = labels;
  
  // Find which indices are to be swapped
  cluster_k = find(loc_labs == kstar);
  cluster_k_prime = find(loc_labs == k_prime);
  
  // Rcpp::Rcout << "\n\nLabels before swapping:\n" << dummy_labels;
  // Rcpp::Rcout << "\nk: " << kstar << "\nk': " << k_prime;
  // Rcpp::Rcout << "\nk cluster indices:\n" << cluster_k.t();
  // Rcpp::Rcout << "\n\nk' cluster indices:\n" << cluster_k_prime.t();
  // Rcpp::Rcout << "\nlabels before update:\n" << loc_labs.t();

  // Swap the label associated with the two clusters
  loc_labs.elem(cluster_k).fill(k_prime);
  loc_labs.elem(cluster_k_prime).fill(kstar);

  
  dummy_labels.col(lstar) = loc_labs;
  // Rcpp::Rcout << "\nlabels after update:\n" << loc_labs.t();
  // Rcpp::Rcout << "\n\nLabels after swapping:\n" << dummy_labels;
  
  return dummy_labels;
}

// Check if labels should be swapped to improve correlation of clustering
// across datasets via random sampling.
void mdi::updateLabels() {
  if(L == 1) {
    return;
  }
  for(uword l = 0; l < L; l++) {
    updateLabelsViewL(l);
  }
};

void mdi::updateLabelsViewL(uword lstar) {
  
  bool multipleUnfixedComponents = (K_unfixed(lstar) > 1),
    accept = false,
    not_pertinent_indices = true;
  
  if(! multipleUnfixedComponents) {
    return;
  }
  
  // The other component considered
  uword k_prime = 0;
  
  // Random uniform number
  double 
    
    // The current likelihood
    current_score = 0.0,
    
    // The competitor
    proposed_score = 0.0,
    
    // The acceptance probability
    acceptance_prob = 1.0,
    log_acceptance_prob = 0.0,
    
    // The weight of the kth component if we do accept the swap
    old_weight = 0.0;
  
  // Membership of the kth and k'th components in the lth dataset
  uvec members_lk, members_lk_prime;
  
  // Vector of entries equal to 1/(K - 1) (as we exclude the current label) and
  // its cumulative sum, used to sample another label to consider swapping.
  vec K_inv, K_inv_cum;
  
  umat swapped_labels(N, L);
  K_inv = ones<vec>(K_unfixed(lstar) - 1) * (1.0 / (double)(K_unfixed(lstar) - 1));
  K_inv_cum = cumsum(K_inv);
  
  for(uword k = K_fixed(lstar); k < K(lstar); k++) {
    
    // The score associated with the current labelling
    current_score = calcScore(lstar, labels);
    
    // Select another label randomly
    k_prime = sampleLabel(k - K_fixed(lstar), K_inv_cum) + K_fixed(lstar);
    if(k_prime > (K(lstar) - 1)) {
      Rcpp::stop("\nk_prime exceeds the number of components.\n");
    }
    if(k_prime == k){
      Rcpp::stop("\nk_prime equals k.\n");
    }
    
    not_pertinent_indices = (N_k(k, lstar) == 0) && (N_k(k_prime, lstar) == 0);
    if(not_pertinent_indices) {
      continue;
    }
    
    // The label matrix updated with the swapped labels
    swapped_labels = swapLabels(lstar, k, k_prime);
    
    // The score for the swap
    proposed_score = calcScore(lstar, swapped_labels);
    
    // The log acceptance probability
    log_acceptance_prob = proposed_score - current_score;

    // Rcpp::Rcout << "\ncurr prob: " << exp(current_score);
    // Rcpp::Rcout << "\nalt prob: " << exp(proposed_score);
    // Rcpp::Rcout << "\nlog acceptance: " << log_acceptance_prob;

    acceptance_prob = std::min(1.0, std::exp(log_acceptance_prob));
    accept = metropolisAcceptanceStep(acceptance_prob);
    
    // If we accept the label swap, update labels, weights and score
    if(accept) {
      acceptance_count++;
      // Rcpp::Rcout << "\n\nLabels before swapping:\n" << labels;
      // Rcpp::Rcout << "\n\nLabels after swapping:\n" << swapped_labels;
      
      // Update the current score
      current_score = proposed_score;
      labels = swapped_labels;
      
      // Pass the new labels from the mixture level back to the MDI level.
      mixtures[lstar]->labels = labels.col(lstar);
      
      // Update the component weights
      old_weight = w(k, lstar);
      w(k, lstar) = w(k_prime, lstar);
      w(k_prime, lstar) = old_weight;

      members_lk = members.slice(lstar).col(k_prime);
      members_lk_prime = members.slice(lstar).col(k);

      members.slice(lstar).col(k) = members_lk;
      members.slice(lstar).col(k_prime) = members_lk_prime;

      N_k(k, lstar) = accu(members_lk);
      N_k(k_prime, lstar) = accu(members_lk_prime);

      // Pass the allocation count down to the mixture
      // (this is used in the parameter sampling)
      mixtures[lstar]->members.col(k) = members_lk;
      mixtures[lstar]->members.col(k_prime) = members_lk_prime;
      mixtures[lstar]->N_k = N_k(span(0, K(lstar) - 1), lstar);
      
    }
  } 
};
