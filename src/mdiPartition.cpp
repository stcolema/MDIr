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

}

double mdiPartitionSum(
    const arma::mat& w,
    const arma::uvec& K,
    const arma::mat& phi
) {
  const uword L = w.n_cols;
  if(L == 0 || L > MAX_VIEWS) {
    Rcpp::stop("mdiPartitionSum: the number of views must be between 1 and %d.", (int) MAX_VIEWS);
  }
  const uword n_sets = (uword) 1 << L, full = n_sets - 1;

  std::vector<double> P(n_sets, 1.0), C(n_sets, 0.0), s(n_sets, 0.0), G(n_sets, 0.0);

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

  // s_X: sum over components present in every view of X of the weight product
  for(uword X = 1; X < n_sets; X++) {
    uword K_min = w.n_rows;
    for(uword l = 0; l < L; l++) {
      if(X & ((uword) 1 << l)) {
        K_min = std::min(K_min, (uword) K(l));
      }
    }
    double total = 0.0;
    for(uword k = 0; k < K_min; k++) {
      double prod_w = 1.0;
      for(uword l = 0; l < L; l++) {
        if(X & ((uword) 1 << l)) {
          prod_w *= w(k, l);
        }
      }
      total += prod_w;
    }
    s[X] = total;
  }

  // C_X (connected spanning edge sets), then G_X (sum over set partitions of X)
  // Proper subsets always have smaller numeric value, so a forward pass works.
  G[0] = 1.0;
  for(uword X = 1; X < n_sets; X++) {
    const uword low = lowBit(X), rest = X ^ low;
    double c = P[X], g = 0.0;

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

    sub = rest;
    while(true) {
      const uword A = sub | low;
      g += C[A] * s[A] * G[X ^ A];
      if(sub == 0) break;
      sub = (sub - 1) & rest;
    }
    G[X] = g;
  }
  return G[full];
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
