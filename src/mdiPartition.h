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
// Independent of the recursion behind mdiWeightRates(), which the sampler uses; kept as the
// reference against which the tests check it.
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

// Log of the Metropolis-Hastings ratio for exchanging components k and kprime of
// view lstar: their labels, weights and component parameters move together. The
// weights and parameters have exchangeable priors within a view and the data
// likelihood and the products of weights over items are carried along by the
// exchange, so only two terms change:
//
//   sum_{m != lstar} log(1 + phi(m, lstar)) (A' - A)  -  v (Z' - Z),
//
// where A (A') counts the items whose label in view m equals their label in view
// lstar before (after) the exchange and Z' is Z with w(k, lstar) and
// w(kprime, lstar) interchanged. Only the weights of view lstar are interchanged:
// Z is unchanged by the same permutation of every view (a relabelling of the
// joint components), so exchanging the other views' weights would drop the
// change in Z that the move produces. Z_current is Z for w; Z' is returned in
// Z_swapped.
double mdiSwapLogRatio(
    const arma::umat& labels,
    const arma::mat& phi,
    const arma::mat& w,
    const arma::uvec& K,
    const std::vector<double>& C,
    double v,
    arma::uword lstar,
    arma::uword k,
    arma::uword kprime,
    double Z_current,
    double& Z_swapped
);

// === Marginal likelihood of an item and its class probabilities ==============
//
// For an item with per-view component likelihoods g(k, l) = p(x_l | component k
// of view l) (outliers marginalised where a view has them), summing the joint
// density over all component assignments gives the same multilinear polynomial
// as Z with the weights w(k, l) replaced by w(k, l) g(k, l):
//
//   p(x | w, phi, theta) = Z(w * g) / Z(w).
//
// Each column of g is rescaled by its maximum before the sum so that nothing
// underflows; the rescaling is added back on the log scale. A column may hold
// -Inf entries (a component an item cannot belong to, e.g. an observed label).

// log Z(w * g), the log of the numerator above. log_g is K_max x L on the log
// scale. Returns -Inf if some view has no component with positive likelihood.
double mdiLogNumerator(
    const arma::mat& w,
    const arma::uvec& K,
    const std::vector<double>& C,
    const arma::mat& log_g
);

// As mdiLogNumerator(), and also the K_max x L matrix whose entry (k, l) is the
// posterior probability that the item belongs to component k of view l given its
// data, p(c_l = k | x) = w(k, l) g(k, l) dZ(w * g) / dw(k, l) / Z(w * g). Entries
// beyond K(l) are zero. The columns sum to one.
arma::mat mdiClassProbabilities(
    const arma::mat& w,
    const arma::uvec& K,
    const std::vector<double>& C,
    const arma::mat& log_g,
    double& log_numerator
);

// === Collapsed conditional of a phi ==========================================
//
// With the strategic latent variable integrated out, the conditional of phi(l, m)
// given the weights and the labels is
//
//   p(phi) propto phi^(shape - 1) exp(-rate phi) (1 + phi)^N_lm (A + B phi)^(-N),
//
// because Z is linear in phi(l, m): Z = A + B phi with A = Z at phi(l, m) = 0
// and B = dZ / dphi(l, m). Draw from it with one slice-sampling update (Neal,
// 2003) on log(phi) from the current value. The update leaves the conditional
// invariant whatever the initial width, so nothing needs tuning; the width only
// affects the cost of the draw, and each evaluation is O(1) once A and B are
// known.
double mdiSamplePhiSlice(
    double phi,
    double N_lm,
    double N,
    double A,
    double B,
    double shape,
    double rate
);

// Log of the unnormalised collapsed conditional in u = log(phi), including the
// Jacobian (used by the sampler above and by tests).
double mdiLogPhiConditional(
    double u,
    double N_lm,
    double N,
    double A,
    double B,
    double shape,
    double rate
);

// === Exact joint draw of the labels of one item across views ================
//
// Given per-view scores a_l(k) (the component weight times the likelihood of the item),
// the joint conditional of the labels (c_1, ..., c_L) of one item is
//
//   p(c) proportional to prod_l a_l(c_l) prod_{l<m} (1 + phi_lm 1[c_l = c_m]).
//
// It is drawn view by view without enumerating the K_1 x ... x K_L combinations. The
// marginal of the first view is a_1(k) dZ(a) / da_1(k) (the rates of mdiWeightRates());
// fixing c_1 = k multiplies a_m(k) by (1 + phi_1m) in every later view m, and the later views
// again have this form with the sub-matrix of phi. The cost is O(L 3^L + L K) per item.

// C tables for every suffix of the views: element t holds mdiConnectedSums() of the
// views t, ..., L - 1. Depends on phi only, so it can be reused for every item.
std::vector< std::vector<double> > mdiSuffixConnectedSums(const arma::mat& phi);

// Draw the labels of the views of one block. G is K_max x b (a_l(k), zero beyond K(l) or
// where the component is excluded), K and phi are restricted to the block, and suffix_C is
// mdiSuffixConnectedSums(phi). Each column of G is rescaled internally.
arma::uvec mdiSampleJointLabels(
    arma::mat G,
    const arma::uvec& K,
    const arma::mat& phi,
    const std::vector< std::vector<double> >& suffix_C
);

// The labels of a block of views for one item, given log_g (K_max x L, log of the likelihood
// of the item in each component, -Inf where excluded), the weights w and the labels
// `current` of the item in every view. Views outside the block keep their labels and
// enter through the factors (1 + phi) they give to the matching component of the block's
// views. Returns the new labels of the block's views, in the order of `block`.
// `suffix_C` may be null, in which case the tables are built here.
arma::uvec mdiSampleJointBlock(
    const arma::mat& log_g,
    const arma::mat& w,
    const arma::uvec& K,
    const arma::mat& phi,
    const arma::uvec& block,
    const arma::uvec& current,
    const std::vector< std::vector<double> >* suffix_C
);

#endif /* MDIPARTITION_H */
