// mdiPartition.cpp
// =============================================================================
# include "mdiPartition.h"

using namespace arma ;

namespace {

// Number of views above which the O(3^L) tables are refused (3^16 ~ 4.3e7).
const uword MAX_VIEWS = 16;

inline uword lowBit(uword x) {
  return x & (~x + 1);
}

// Index of the highest set bit of x (x > 0)
inline uword highIndex(uword x) {
  uword h = 0;
  while((x >> (h + 1)) != 0) {
    h++;
  }
  return h;
}

void checkNumberOfViews(uword L) {
  if(L == 0 || L > MAX_VIEWS) {
    Rcpp::stop("mdiPartitionSum: the number of views must be between 1 and %d.", (int) MAX_VIEWS);
  }
}

// Fills K_min[X], the smallest number of components over the views in X
// (w.n_rows for the empty set), and the index of the highest view in X.
void subsetTables(
    const arma::mat& w,
    const arma::uvec& K,
    std::vector<uword>& K_min,
    std::vector<uword>& high
) {
  const uword L = w.n_cols, n_sets = (uword) 1 << L;
  K_min.assign(n_sets, (uword) w.n_rows);
  high.assign(n_sets, 0);
  for(uword X = 1; X < n_sets; X++) {
    const uword h = highIndex(X);
    high[X] = h;
    K_min[X] = std::min(K_min[X ^ ((uword) 1 << h)], (uword) K(h));
  }
}

// s_X: sum over the components present in every view of X of the product of
// their weights. The products are built by adding the views of X in increasing
// order, the same order as a direct loop over the views.
void weightProductSums(
    const arma::mat& w,
    const std::vector<uword>& K_min,
    const std::vector<uword>& high,
    std::vector<double>& s
) {
  const uword L = w.n_cols, n_sets = (uword) 1 << L;
  s.assign(n_sets, 0.0);
  std::vector<double> prod_w(n_sets, 1.0);
  for(uword k = 0; k < w.n_rows; k++) {
    for(uword X = 1; X < n_sets; X++) {
      const uword h = high[X];
      prod_w[X] = prod_w[X ^ ((uword) 1 << h)] * w(k, h);
      if(k < K_min[X]) {
        s[X] += prod_w[X];
      }
    }
  }
}

// G_X: sum over the set partitions of X of prod_{B} C_B s_B. G[full] is Z.
void setPartitionSums(
    const std::vector<double>& C,
    const std::vector<double>& s,
    std::vector<double>& G
) {
  const uword n_sets = C.size();
  std::vector<double> Cs(n_sets);
  for(uword X = 0; X < n_sets; X++) {
    Cs[X] = C[X] * s[X];
  }
  G.assign(n_sets, 0.0);
  G[0] = 1.0;
  for(uword X = 1; X < n_sets; X++) {
    const uword low = lowBit(X), rest = X ^ low;
    double g = 0.0;
    uword sub = rest;
    while(true) {
      const uword A = sub | low;
      g += Cs[A] * G[X ^ A];
      if(sub == 0) break;
      sub = (sub - 1) & rest;
    }
    G[X] = g;
  }
}

}

std::vector<double> mdiConnectedSums(const arma::mat& phi) {
  const uword L = phi.n_cols;
  checkNumberOfViews(L);
  const uword n_sets = (uword) 1 << L;

  std::vector<double> P(n_sets, 1.0), C(n_sets, 0.0);

  // P_X: product of (1 + phi) over the pairs within X
  for(uword X = 1; X < n_sets; X++) {
    const uword low = lowBit(X), rest = X ^ low;
    uword l = 0;
    while(((uword) 1 << l) != low) { l++; }
    double p = P[rest];
    for(uword m = l + 1; m < L; m++) {
      if(rest & ((uword) 1 << m)) {
        p *= (1.0 + phi(l, m));
      }
    }
    P[X] = p;
  }

  // C_X (connected spanning edge sets). Proper subsets always have smaller
  // numeric value, so a forward pass works.
  for(uword X = 1; X < n_sets; X++) {
    const uword low = lowBit(X), rest = X ^ low;
    double c = P[X];

    // Enumerate A = low | sub for every sub of rest
    uword sub = rest;
    while(true) {
      const uword A = sub | low;
      if(A != X) {
        c -= C[A] * P[X ^ A];
      }
      if(sub == 0) break;
      sub = (sub - 1) & rest;
    }
    C[X] = c;
  }
  return C;
}

double mdiPartitionSumFromC(
    const arma::mat& w,
    const arma::uvec& K,
    const std::vector<double>& C
) {
  const uword L = w.n_cols;
  checkNumberOfViews(L);
  if(C.size() != ((uword) 1 << L)) {
    Rcpp::stop("mdiPartitionSum: the connected sums do not match the number of views.");
  }
  std::vector<uword> K_min, high;
  std::vector<double> s, G;
  subsetTables(w, K, K_min, high);
  weightProductSums(w, K_min, high, s);
  setPartitionSums(C, s, G);
  return G[((uword) 1 << L) - 1];
}

double mdiPartitionSum(
    const arma::mat& w,
    const arma::uvec& K,
    const arma::mat& phi
) {
  checkNumberOfViews(w.n_cols);
  return mdiPartitionSumFromC(w, K, mdiConnectedSums(phi));
}

arma::vec mdiWeightRates(
    const arma::mat& w,
    const arma::uvec& K,
    const std::vector<double>& C,
    arma::uword lstar
) {
  const uword L = w.n_cols;
  checkNumberOfViews(L);
  if(C.size() != ((uword) 1 << L)) {
    Rcpp::stop("mdiWeightRates: the connected sums do not match the number of views.");
  }
  const uword n_sets = (uword) 1 << L, full = n_sets - 1, lbit = (uword) 1 << lstar;
  const uword rest_mask = full ^ lbit;

  std::vector<uword> K_min, high;
  std::vector<double> s, G;
  subsetTables(w, K, K_min, high);
  weightProductSums(w, K_min, high, s);
  setPartitionSums(C, s, G);

  // Z is linear in w(kstar, lstar). Collecting the terms by the block B of the
  // set partition that holds lstar gives
  //   dZ / dw(k, lstar) = sum_B C_B prod_{j in B \ lstar} w(k, j) G[full \ B],
  // over the blocks B in which every other view has more than k components.
  // G[full \ B] does not involve view lstar.
  arma::vec rates(K(lstar), arma::fill::zeros);
  std::vector<double> prod_w(n_sets, 1.0);
  for(uword k = 0; k < K(lstar); k++) {
    for(uword X = 1; X < n_sets; X++) {
      if(X & lbit) {
        continue;
      }
      const uword h = high[X];
      prod_w[X] = prod_w[X ^ ((uword) 1 << h)] * w(k, h);
    }
    double total = 0.0;
    uword sub = rest_mask;
    while(true) {
      if(k < K_min[sub]) {
        const uword B = sub | lbit;
        total += C[B] * prod_w[sub] * G[full ^ B];
      }
      if(sub == 0) break;
      sub = (sub - 1) & rest_mask;
    }
    rates(k) = total;
  }
  return rates;
}

double mdiWeightRate(
    const arma::mat& w,
    const arma::uvec& K,
    const arma::mat& phi,
    arma::uword lstar,
    arma::uword kstar
) {
  arma::mat w_unit = w;
  w_unit.col(lstar).zeros();
  w_unit(kstar, lstar) = 1.0;
  return mdiPartitionSum(w_unit, K, phi);
}

double mdiPhiRate(
    const arma::mat& w,
    const arma::uvec& K,
    const arma::mat& phi,
    arma::uword l,
    arma::uword m
) {
  const uword L = w.n_cols;

  // Merged problem: the merged view sits in position 0
  arma::mat w_merged(w.n_rows, L - 1, arma::fill::zeros);
  arma::mat phi_merged(L - 1, L - 1, arma::fill::zeros);
  arma::uvec K_merged(L - 1);

  const uword K_lm = std::min((uword) K(l), (uword) K(m));
  for(uword k = 0; k < K_lm; k++) {
    w_merged(k, 0) = w(k, l) * w(k, m);
  }
  K_merged(0) = K_lm;

  std::vector<uword> others;
  for(uword j = 0; j < L; j++) {
    if(j != l && j != m) {
      others.push_back(j);
    }
  }
  for(uword a = 0; a < others.size(); a++) {
    const uword j = others[a];
    w_merged.col(a + 1) = w.col(j);
    K_merged(a + 1) = K(j);
    const double up = (1.0 + phi(l, j)) * (1.0 + phi(m, j)) - 1.0;
    phi_merged(0, a + 1) = up;
    phi_merged(a + 1, 0) = up;
    for(uword b = a + 1; b < others.size(); b++) {
      phi_merged(a + 1, b + 1) = phi(j, others[b]);
      phi_merged(b + 1, a + 1) = phi(j, others[b]);
    }
  }
  return mdiPartitionSum(w_merged, K_merged, phi_merged);
}
