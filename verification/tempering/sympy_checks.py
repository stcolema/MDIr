"""Symbolic and exact-arithmetic checks for the tempered sampler and replica exchange.

Run:  python3 sympy_checks.py
Requires sympy and mpmath. Every check prints PASS or raises.

1. Tempered normal-inverse-gamma update (scalar) and its marginal likelihood.
2. Tempered normal-inverse-Wishart update (P = 2, symbolic entries).
3. Tempered Dirichlet-categorical update and marginal likelihood.
4. The exchange acceptance ratio: the prior cancels.
5. Exact rational check on a small finite product space that
   - each swap kernel, each even/odd round, the DEO and SEO iterations leave the
     product target invariant, and
   - the DEO two-round kernel is NOT reversible (invariant but non-reversible),
   - a local kernel product followed by exchanges is invariant.
"""
import itertools
import random

import mpmath as mp
import sympy as sp

mp.mp.dps = 30


def ok(name):
    print(f"PASS  {name}")


# ---------------------------------------------------------------- 1. scalar NIG
def check_nig():
    mu, s2 = sp.symbols("mu sigma2", positive=True, real=True)
    beta, kappa, nu, psi, xi = sp.symbols("beta kappa nu psi xi", positive=True, real=True)
    n, S1, S2 = sp.symbols("n S1 S2", positive=True, real=True)
    xbar = S1 / n
    SS = S2 - S1**2 / n  # centred sum of squares

    # log of the tempered likelihood, prior, as functions of (mu, sigma2).
    # likelihood^beta: prod_i N(x_i | mu, s2)^beta, using sum x = S1, sum x^2 = S2
    sum_sq = S2 - 2 * mu * S1 + n * mu**2
    log_lik_b = beta * (-n / 2 * sp.log(2 * sp.pi * s2) - sum_sq / (2 * s2))
    log_prior_mu = -sp.log(2 * sp.pi * s2 / kappa) / 2 - kappa * (mu - xi) ** 2 / (2 * s2)
    log_prior_s2 = (nu / 2) * sp.log(psi / 2) - sp.loggamma(nu / 2) - (nu / 2 + 1) * sp.log(s2) - psi / (2 * s2)
    lhs = log_lik_b + log_prior_mu + log_prior_s2

    kn = kappa + beta * n
    nun = nu + beta * n
    mun = (kappa * xi + beta * S1) / kn
    psin = psi + beta * SS + kappa * beta * n / kn * (xbar - xi) ** 2
    # claimed posterior: N(mu | mun, s2/kn) IG(s2 | nun/2, psin/2) times marginal likelihood m_beta
    log_post_mu = -sp.log(2 * sp.pi * s2 / kn) / 2 - kn * (mu - mun) ** 2 / (2 * s2)
    log_post_s2 = (nun / 2) * sp.log(psin / 2) - sp.loggamma(nun / 2) - (nun / 2 + 1) * sp.log(s2) - psin / (2 * s2)
    log_m = (
        -beta * n / 2 * sp.log(2 * sp.pi)
        + sp.log(kappa / kn) / 2
        + (nu / 2) * sp.log(psi / 2)
        - sp.loggamma(nu / 2)
        + sp.loggamma(nun / 2)
        - (nun / 2) * sp.log(psin / 2)
    )
    diff = sp.simplify(sp.expand_log(lhs - (log_post_mu + log_post_s2 + log_m), force=True))
    assert sp.simplify(diff) == 0, diff
    ok("tempered NIG: prior x L^beta = posterior x m_beta, symbolically (identity in mu, sigma2)")

    # same, with the claimed updates being exactly the untempered update at n_eff = beta n, scatter beta*SS
    n_eff = beta * n
    assert sp.simplify(kn - (kappa + n_eff)) == 0
    assert sp.simplify(nun - (nu + n_eff)) == 0
    assert sp.simplify(kappa * n_eff / (kappa + n_eff) - kappa * beta * n / kn) == 0
    ok("tempered NIG update = untempered update with n_eff = beta n and scatter beta*SS")

    # numeric check of the marginal likelihood by 2-d quadrature
    vals = {beta: 0.37, kappa: 0.8, nu: 3.0, psi: 1.7, xi: 0.4, n: 6, S1: 3.1, S2: 9.9}
    f = sp.lambdify((mu, s2), sp.exp(lhs.subs(vals)), "mpmath")
    num = mp.quad(lambda m_: mp.quad(lambda v_: f(m_, v_), [0, 0.5, 2, 10, mp.inf]), [-mp.inf, -2, 0, 2, mp.inf])
    ana = mp.e ** sp.N(log_m.subs(vals), 30)
    assert abs(num / ana - 1) < 1e-8, (num, ana)
    ok(f"marginal likelihood m_beta matches 2-d quadrature (rel. err {float(abs(num / ana - 1)):.1e})")


# ---------------------------------------------------------------- 2. NIW, P = 2
def check_niw():
    P = 2
    beta, kappa, nu, n = sp.symbols("beta kappa nu n", positive=True, real=True)
    # symbols for the sufficient statistics: sum x (vector), sum x x^T (symmetric)
    s1 = sp.Matrix(sp.symbols("s1_0 s1_1", real=True))
    a, b, c = sp.symbols("a b c", real=True)
    S2 = sp.Matrix([[a, b], [b, c]])
    xi = sp.Matrix(sp.symbols("xi0 xi1", real=True))
    m = sp.Matrix(sp.symbols("m0 m1", real=True))
    Sig = sp.Matrix([[sp.Symbol("v11", positive=True), sp.Symbol("v12", real=True)],
                     [sp.Symbol("v12", real=True), sp.Symbol("v22", positive=True)]])
    Psi = sp.Matrix([[sp.Symbol("p11", positive=True), sp.Symbol("p12", real=True)],
                     [sp.Symbol("p12", real=True), sp.Symbol("p22", positive=True)]])
    Si = Sig.inv()
    # exponent (the part of the log density that depends on mu, up to terms in Sigma only)
    # beta * sum_i (x_i - m)^T Si (x_i - m) + kappa (m - xi)^T Si (m - xi)
    tr_term = (Si * S2).trace()
    lik_quad = tr_term - 2 * (m.T * Si * s1)[0] + n * (m.T * Si * m)[0]
    lhs = beta * lik_quad + kappa * ((m - xi).T * Si * (m - xi))[0]

    kn = kappa + beta * n
    mn = (kappa * xi + beta * s1) / kn
    xbar = s1 / n
    scatter = S2 - s1 * s1.T / n
    # claimed: (m - mn)^T Si (m - mn) kn + tr(Si * (beta*scatter + kappa beta n / kn (xbar - xi)(xbar - xi)^T))
    rhs = kn * ((m - mn).T * Si * (m - mn))[0] + (
        Si * (beta * scatter + (kappa * beta * n / kn) * (xbar - xi) * (xbar - xi).T)
    ).trace()
    d = sp.simplify(sp.expand(lhs - rhs))
    assert d == 0, d
    ok("tempered NIW (P = 2): completing the square, symbolic identity in (mu, Sigma)")
    # The log|Sigma| terms. With ld = log|Sigma| and P = 2:
    #   prior IW(Psi, nu):            -(nu + P + 1)/2 * ld
    #   prior N(mu | xi, Sigma/kappa): -1/2 * ld  (+ constants in kappa)
    #   likelihood^beta:               -beta n / 2 * ld
    # claimed posterior IW(Psi_n, nu + beta n) x N(mu | mn, Sigma / kn): -(nu + beta n + P + 1)/2 * ld - 1/2 * ld
    ld = sp.Symbol("ld")
    lhs_ld = -(nu + P + 1) / 2 * ld - sp.Rational(1, 2) * ld - beta * n / 2 * ld
    rhs_ld = -(nu + beta * n + P + 1) / 2 * ld - sp.Rational(1, 2) * ld
    assert sp.simplify(lhs_ld - rhs_ld) == 0
    # constants: the normaliser ratio (kappa/kn)^(P/2) is free of (mu, Sigma), so it
    # only enters the marginal likelihood
    ok("tempered NIW degrees of freedom nu + beta n: the log|Sigma| coefficients agree (symbolic)")


# ---------------------------------------------------------------- 3. Dirichlet
def check_dirichlet():
    alpha = [0.7, 1.3, 0.5]
    counts = [3, 1, 4]
    beta = 0.43

    def logB(v):
        return sum(mp.loggamma(x) for x in v) - mp.loggamma(sum(v))

    # integrate over the 2-simplex
    def integrand(t1, t2):
        t3 = 1 - t1 - t2
        if t3 <= 0:
            return mp.mpf(0)
        th = [t1, t2, t3]
        val = mp.mpf(1)
        for ti, ai, ci in zip(th, alpha, counts):
            val *= ti ** (ai - 1 + beta * ci)
        return val

    num = mp.quad(lambda t1: mp.quad(lambda t2: integrand(t1, t2), [0, 1 - t1]), [0, 1])
    ana = mp.e ** (logB([a + beta * c for a, c in zip(alpha, counts)]))
    assert abs(num / ana - 1) < 1e-6, (num, ana)
    ok("tempered Dirichlet: integral of prod theta^(alpha - 1 + beta n) = B(alpha + beta n) (quadrature)")


# ---------------------------------------------------------------- 4. exchange ratio
def check_exchange_ratio():
    bi, bj, li, lj, pi_, pj = sp.symbols("beta_i beta_j ell_i ell_j p_i p_j", real=True)
    # state i at temperature i has log-lik li and prior weight pi_, similarly j
    num = sp.exp(bi * lj) * pj * sp.exp(bj * li) * pi_  # states exchanged
    den = sp.exp(bi * li) * pi_ * sp.exp(bj * lj) * pj
    r = sp.simplify(num / den)
    assert sp.simplify(r - sp.exp((bi - bj) * (lj - li))) == 0
    ok("exchange ratio = exp((beta_i - beta_j)(ell_j - ell_i)); priors cancel")


# ---------------------------------------------------------------- 5. exact finite chain
def check_exact_pt():
    random.seed(1)
    S = [0, 1, 2]
    T = 3
    # integer 'inverse temperatures' k_t and rational base weights keep every
    # quantity exactly rational: pi_t(s) = p(s) * q(s)^k_t, with q the likelihood
    p = {0: sp.Rational(1, 2), 1: sp.Rational(3, 10), 2: sp.Rational(1, 5)}
    q = {0: sp.Rational(1, 7), 1: sp.Rational(3, 1), 2: sp.Rational(1, 2)}
    ks = [0, 1, 2]

    def pit(t, s):
        return p[s] * q[s] ** ks[t]

    states = list(itertools.product(S, repeat=T))
    idx = {s: i for i, s in enumerate(states)}
    Nn = len(states)
    Pi = sp.Matrix([sp.prod([pit(t, s[t]) for t in range(T)]) for s in states])

    def invariant(K, name):
        res = (Pi.T * K - Pi.T)
        assert all(sp.simplify(x) == 0 for x in res), name
        # rows sum to one
        assert all(sp.simplify(sum(K.row(i)) - 1) == 0 for i in range(Nn)), name + " (rows)"

    def reversible(K):
        for i in range(Nn):
            for j in range(Nn):
                if sp.simplify(Pi[i] * K[i, j] - Pi[j] * K[j, i]) != 0:
                    return False
        return True

    def swap_kernel(i, j):
        K = sp.zeros(Nn, Nn)
        for s in states:
            s2 = list(s)
            s2[i], s2[j] = s2[j], s2[i]
            s2 = tuple(s2)
            ratio = Pi[idx[s2]] / Pi[idx[s]]
            a = sp.Min(1, ratio)
            K[idx[s], idx[s2]] += a
            K[idx[s], idx[s]] += 1 - a
        return K

    # local kernel at temperature t: Metropolis on the three states with a
    # cyclic proposal (not reversible on its own) plus the reverse, i.e. a
    # symmetric random-walk Metropolis on the 3-cycle, exactly rational
    def local_kernel_state(t):
        L = sp.zeros(3, 3)
        for s in S:
            for d in (1, -1):
                s2 = (s + d) % 3
                a = sp.Min(1, pit(t, s2) / pit(t, s))
                L[s, s2] += sp.Rational(1, 2) * a
                L[s, s] += sp.Rational(1, 2) * (1 - a)
        return L

    Ls = [local_kernel_state(t) for t in range(T)]
    # each local kernel is pi_t-invariant (exact)
    for t in range(T):
        v = sp.Matrix([pit(t, s) for s in S]).T * Ls[t] - sp.Matrix([pit(t, s) for s in S]).T
        assert all(sp.simplify(x) == 0 for x in v)
    ok("local Metropolis kernels are invariant for their tempered targets (exact)")

    Kloc = sp.zeros(Nn, Nn)
    for s in states:
        for s2 in states:
            Kloc[idx[s], idx[s2]] = sp.prod([Ls[t][s[t], s2[t]] for t in range(T)])
    invariant(Kloc, "product of local kernels")
    ok("product of local kernels is invariant for the product target (exact)")

    K01, K12 = swap_kernel(0, 1), swap_kernel(1, 2)
    for name, K in (("swap(0,1)", K01), ("swap(1,2)", K12)):
        invariant(K, name)
        assert reversible(K), name
    ok("each swap kernel is invariant and reversible (exact)")

    Keven, Kodd = K01, K12      # with T = 3 the even pairs are {(0,1)} and the odd {(1,2)}
    KDEO = Keven * Kodd          # two consecutive DEO rounds (even, then odd)
    KSEO = sp.Rational(1, 2) * Keven + sp.Rational(1, 2) * Kodd
    invariant(KDEO, "DEO two-round")
    invariant(KSEO, "SEO")
    ok("DEO (even then odd) and SEO (random parity) exchange kernels are invariant (exact)")
    assert not reversible(KDEO)
    assert reversible(KSEO)
    ok("DEO two-round kernel is invariant but NOT reversible; SEO is reversible (exact)")

    Kfull = Kloc * Keven * Kloc * Kodd
    invariant(Kfull, "local + DEO iteration")
    ok("local sweep + even round + local sweep + odd round is invariant (exact)")

    # marginal at the coldest temperature is pi_cold (normalised)
    z = [sum(pit(t, s) for s in S) for t in range(T)]
    marg = [sum(Pi[idx[s]] for s in states if s[T - 1] == v) / sp.prod(z) for v in S]
    assert all(sp.simplify(marg[v] - pit(T - 1, v) / z[T - 1]) == 0 for v in S)
    ok("cold-coordinate marginal of the product target is the normalised pi_cold (exact)")


if __name__ == "__main__":
    check_nig()
    check_niw()
    check_dirichlet()
    check_exchange_ratio()
    check_exact_pt()
    print("all symbolic / exact checks passed")
