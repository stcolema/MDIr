// mdiPartition.h
// =============================================================================
// Exact evaluation of the MDI normalising constant and the quantities derived
// from it (rates of the conditional distributions of the component weights and
// the phi parameters) without enumerating all K^L joint component assignments.
//
// For views l = 1, ..., L with component weights w(k, l) and pairwise
// parameters phi(l, m) the MDI normalising constant is
//
//   Z = sum_{k_1..k_L} prod_l w(k_l, l) prod_{l<m} (1 + phi(l, m) 1[k_l = k_m]).
//
// Expanding the product over pairs into a sum over edge subsets S, each S
// forces the views in a connected component of S to share a component index, so
//
//   Z = sum over set partitions pi of {1..L} of prod_{B in pi} C_B s_B,
//
// where s_B = sum_k prod_{l in B} w(k, l) (k restricted to components present
// in every view of B) and C_B is the sum over connected spanning edge sets of B
// of prod phi. C_B is obtained from P_B = prod_{l<m in B} (1 + phi(l, m)) via
//
//   P_B = sum_{A subset B, min(B) in A} C_A P_{B \ A},
//
// and the sum over set partitions with the same recursion. The cost is
// O(3^L + 2^L K) rather than O(K^L L^2).
//
// C_B depends on the phis alone, so it is computed once per set of phis
// (mdiConnectedSums) and reused for every evaluation that shares them: the
// normalising constant, the rates of all the weights of a view and the label
// swaps.
#ifndef MDIPARTITION_H
#define MDIPARTITION_H

# include <RcppArmadillo.h>
# include <vector>

// C_X for every subset X of the views (bit l of X set when view l is in X),
// from the symmetric L x L matrix of phis (the diagonal is ignored). Depends on
// phi only.
std::vector<double> mdiConnectedSums(const arma::mat& phi);

// Z from precomputed connected sums (see mdiConnectedSums).
double mdiPartitionSumFromC(
    const arma::mat& w,
    const arma::uvec& K,
    const std::vector<double>& C
);

// dZ / dw(k, lstar) for every k < K(lstar) from one pass over the set
// partitions. Equal, up to rounding, to mdiWeightRate() for each k.
arma::vec mdiWeightRates(
    const arma::mat& w,
    const arma::uvec& K,
    const std::vector<double>& C,
    arma::uword lstar
);

// Z for weights w (K_max x L), numbers of components K (length L) and the
// symmetric L x L matrix of phis (the diagonal is ignored).
double mdiPartitionSum(
    const arma::mat& w,
    const arma::uvec& K,
    const arma::mat& phi
);

// dZ / dw(kstar, lstar) for a single weight: the rate of the conditional for w(kstar, lstar) is
// w_rate_prior + v * this. Z is multilinear in the columns of w, so this is Z
// evaluated with column lstar replaced by the kstar-th unit vector.
double mdiWeightRate(
    const arma::mat& w,
    const arma::uvec& K,
    const arma::mat& phi,
    arma::uword lstar,
    arma::uword kstar
);

// dZ / dphi(l, m): Z is linear in phi(l, m), and its coefficient is the
// normalising constant of the problem in which views l and m are merged into
// one view with weights w(k, l) w(k, m) and phi'(j) = (1 + phi(l, j))
// (1 + phi(m, j)) - 1 against every other view j.
double mdiPhiRate(
    const arma::mat& w,
    const arma::uvec& K,
    const arma::mat& phi,
    arma::uword l,
    arma::uword m
);

#endif /* MDIPARTITION_H */
