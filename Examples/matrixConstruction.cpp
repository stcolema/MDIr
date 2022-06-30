
# include <RcppArmadillo.h>

using namespace arma ;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat firstCovProduct(uword N, uword P, mat entries) {
  uword jj_bound = 0;
  double new_entry = 0.0, lambda = 0.0;
  uvec 
    N_inds(N),
    P_inds(P), 
    current_elements(4),
    first_element(1),
    second_element(1),
    third_element(1),
    fourth_element(1),
    j_inds,
    rel_cols_inds(N),
    first_col_inds(N),
    second_col_inds(N),
    third_col_inds(N),
    fourth_col_inds(N);
  
  mat new_mat(P, N * P);
  new_mat.zeros();
  
  current_elements.zeros();

  N_inds = regspace< uvec >(0, N - 1);
  P_inds = regspace< uvec >(0, P - 1);
  
  for(uword ii = 0; ii < P; ii++) {
    
    first_element.fill(ii);
    fourth_element.fill(P - ii - 1);
    
    second_col_inds = regspace< uvec >(ii, P, N * P - 1);
    third_col_inds = regspace< uvec >(P - ii - 1, P, N * P - 1);
    
    jj_bound = std::min(ii + 1, P - ii);
    for(uword jj = 0; jj < jj_bound; jj++) {
      
      second_element.fill(jj);
      third_element.fill(P - jj - 1);
      
      rel_cols_inds = regspace< uvec >(jj, P, P * N - 1);
      first_col_inds = rel_cols_inds;
      fourth_col_inds = regspace< uvec >(P - jj - 1, P, N * P - 1);
      
      new_entry = entries(ii, jj);
      
      new_mat.submat(first_element, first_col_inds).fill(new_entry);
      if(ii != jj) {
        new_mat.submat(second_element, second_col_inds).fill(new_entry);
      }
      if((ii + jj) != (P - 1)) {
        new_mat.submat(third_element, third_col_inds).fill(new_entry);
        new_mat.submat(fourth_element, fourth_col_inds).fill(new_entry);
      }
    }
  }
  return new_mat;
};