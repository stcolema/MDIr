"""Exact (rational arithmetic) enumeration of a tiny SMC sampler.

State space {0, 1, 2}; targets pi_k(s) = p(s) q(s)^k for k = 0..3 (integer powers keep
everything rational; this is the tempering path with integer "inverse temperatures");
kernel at k: random-walk Metropolis on the 3-cycle, invariant for pi_k. N particles.

For every outcome of the algorithm (initial draws, resampling, moves) we sum
probability x estimator, so the expectation is exact. Checks:

 * fixed schedule, any resampling scheme (never / adaptive ESS rule / always;
   multinomial or systematic): E[Zhat] = Z_3 / Z_0 and E[Zhat * fhat(s)] = Z_3 pi_3(s) / Z_0
   for every state s, i.e. the unnormalised estimator is unbiased (Del Moral, 2004, Thm 7.4.2);
 * the self-normalised estimator fhat is NOT unbiased at finite N;
 * an adaptive schedule chosen from the particles' conditional ESS: E[Zhat] != Z (bias),
   which is why exact unbiasedness is claimed for fixed schedules only.

Run: python3 sympy_smc_exact.py
"""
import itertools
from fractions import Fraction as F

S = [0, 1, 2]
p = {0: F(1, 2), 1: F(3, 10), 2: F(1, 5)}
q = {0: F(1, 7), 1: F(3, 1), 2: F(1, 2)}
KMAX = 3


def pi_un(k, s):
    return p[s] * q[s] ** k


Z = [sum(pi_un(k, s) for s in S) for k in range(KMAX + 1)]


def kernel(k):
    M = {}
    for s in S:
        row = {t: F(0) for t in S}
        for d in (1, -1):
            t = (s + d) % 3
            a = min(F(1), pi_un(k, t) / pi_un(k, s))
            row[t] += F(1, 2) * a
            row[s] += F(1, 2) * (1 - a)
        M[s] = row
    return M


M = {k: kernel(k) for k in range(KMAX + 1)}
# invariance of every kernel (exact)
for k in range(KMAX + 1):
    for t in S:
        assert sum(pi_un(k, s) * M[k][s][t] for s in S) == pi_un(k, t)


def resample_outcomes(W, scheme):
    """List of (probability, ancestor tuple) for N = len(W)."""
    N = len(W)
    if scheme == "multinomial":
        out = []
        for anc in itertools.product(range(N), repeat=N):
            pr = F(1)
            for a in anc:
                pr *= W[a]
            if pr > 0:
                out.append((pr, anc))
        return out
    # systematic: one uniform u in [0, 1/N), positions u + j/N
    cum = []
    c = F(0)
    for w in W:
        c += w
        cum.append(c)
    # breakpoints of u in [0, 1/N): where u + j/N crosses a cumulative weight
    bps = {F(0), F(1, N)}
    for ci in cum:
        for j in range(N):
            b = ci - F(j, N)
            if 0 < b < F(1, N):
                bps.add(b)
    bps = sorted(bps)
    out = []
    for a, b in zip(bps[:-1], bps[1:]):
        u = (a + b) / 2
        anc = []
        for j in range(N):
            pos = u + F(j, N)
            i = next(i for i, ci in enumerate(cum) if pos <= ci)
            anc.append(i)
        out.append(((b - a) * N, tuple(anc)))
    return out


def move_outcomes(xs, k):
    """All outcomes of moving every particle independently with kernel k."""
    out = []
    rows = [M[k][x] for x in xs]
    for ys in itertools.product(S, repeat=len(xs)):
        pr = F(1)
        for r, y in zip(rows, ys):
            pr *= r[y]
        if pr > 0:
            out.append((pr, ys))
    return out


def cess(W, ws):
    N = len(W)
    num = sum(Wi * wi for Wi, wi in zip(W, ws)) ** 2
    den = sum(Wi * wi * wi for Wi, wi in zip(W, ws))
    return F(N) * num / den


def expectation(N, schedule_rule, resample_thr, scheme):
    """Return (E[Zhat], {s: E[Zhat * fhat(s)]}, {s: E[fhat(s)]}) exactly.

    schedule_rule: 'fixed' (k = 1, 2, 3) or ('adaptive', alpha): from k jump to k + 2
    when the conditional ESS of that jump is at least alpha * N, else to k + 1.
    resample_thr: resample when ESS < thr * N (0: never; 1: whenever ESS < N).

    Dynamic programme over distinct configurations (k, states, weights). For each we
    carry prob = P(configuration) and pz = E[Zhat so far * 1(configuration)]; both
    updates are linear, so merging configurations is exact.
    """
    levels = {k: {} for k in range(KMAX + 1)}
    pi0 = {s: pi_un(0, s) / Z[0] for s in S}
    for xs in itertools.product(S, repeat=N):
        pr = F(1)
        for x in xs:
            pr *= pi0[x]
        key = (xs, tuple([F(1, N)] * N))
        a, b = levels[0].get(key, (F(0), F(0)))
        levels[0][key] = (a + pr, b + pr)
    for k in range(KMAX):
        for (xs, W), (prob, pz) in levels[k].items():
            if schedule_rule == "fixed":
                k2 = k + 1
            else:
                alpha = schedule_rule[1]
                k2 = k + 1
                if k + 2 <= KMAX:
                    ws2 = [q[x] ** 2 for x in xs]
                    if cess(list(W), ws2) >= alpha * N:
                        k2 = k + 2
            ws = [q[x] ** (k2 - k) for x in xs]
            zinc = sum(Wi * wi for Wi, wi in zip(W, ws))
            Wn = [Wi * wi / zinc for Wi, wi in zip(W, ws)]
            ess = 1 / sum(w * w for w in Wn)
            if ess < resample_thr * N:
                branches = [(pr, [xs[a] for a in anc], [F(1, N)] * N) for pr, anc in resample_outcomes(Wn, scheme)]
            else:
                branches = [(F(1), list(xs), Wn)]
            for pr, xr, Wr in branches:
                for pm, ys in move_outcomes(xr, k2):
                    key = (ys, tuple(Wr))
                    a, b = levels[k2].get(key, (F(0), F(0)))
                    levels[k2][key] = (a + prob * pr * pm, b + pz * zinc * pr * pm)
    ez = F(0)
    ezf = {s: F(0) for s in S}
    ef = {s: F(0) for s in S}
    total = F(0)
    for (xs, W), (prob, pz) in levels[KMAX].items():
        total += prob
        ez += pz
        for s in S:
            fh = sum(Wi for Wi, x in zip(W, xs) if x == s)
            ezf[s] += pz * fh
            ef[s] += prob * fh
    assert total == 1
    return ez, ezf, ef


def check(name, N, rule, thr, scheme, expect_unbiased):
    ez, ezf, ef = expectation(N, rule, thr, scheme)
    target_z = Z[KMAX] / Z[0]
    bias_z = ez - target_z
    biases = {s: ezf[s] - pi_un(KMAX, s) / Z[0] for s in S}
    unb = (bias_z == 0) and all(b == 0 for b in biases.values())
    status = "unbiased" if unb else "biased"
    assert unb == expect_unbiased, (name, bias_z, biases)
    ratio_bias = {s: float(ef[s] - pi_un(KMAX, s) / Z[KMAX]) for s in S}
    print(f"PASS  {name:58s} {status}; E[Zhat] - Z = {float(bias_z):+.3e}; "
          f"self-normalised bias at N={N}: {max(abs(v) for v in ratio_bias.values()):.2e}")


if __name__ == "__main__":
    for scheme in ("multinomial", "systematic"):
        check(f"N=2 fixed schedule, {scheme}, resample always", 2, "fixed", F(2), scheme, True)
        check(f"N=2 fixed schedule, {scheme}, adaptive resampling (ESS < N/2)", 2, "fixed", F(1, 2), scheme, True)
    check("N=2 fixed schedule, never resample (AIS)", 2, "fixed", F(0), "multinomial", True)
    check("N=3 fixed schedule, multinomial, adaptive resampling (ESS < 0.7 N)", 3, "fixed", F(7, 10), "multinomial", True)
    check("N=3 fixed schedule, systematic, adaptive resampling (ESS < 0.7 N)", 3, "fixed", F(7, 10), "systematic", True)
    for alpha in (F(7, 10), F(9, 10)):
        check(f"N=2 ADAPTIVE schedule (conditional ESS >= {float(alpha)} N), multinomial", 2, ("adaptive", alpha), F(1, 2),
              "multinomial", False)
    print("all exact SMC checks passed")
