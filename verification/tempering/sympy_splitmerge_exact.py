"""Exact rational check of the sequentially-allocated split-merge kernel (Dahl 2005 style).

Items 0..n-1, two anchors a, b (items 0, 1) in components A, B; the other items are visited
in a fixed order (the shuffle is part of the proposal; we take one order, then average over
all orders, both exact). Target: pi(c) proportional to prod_n w[c_n] * prod_k m(block_k),
with m(block) = a rational "marginal" that depends only on the block (here a Dirichlet-
multinomial-like product of rising factorials, so no floating point enters).

Proposal: sequentially allocate each free item to A or B with probability proportional to
w[c] * m(block + item)/m(block). The proposal density of a path is prod of those
probabilities; the acceptance is min(1, [pi(new) q(old|new-state)] / [pi(old) q(new|old-state)]).
We verify, in Fraction arithmetic:
  (1) q is a probability distribution over the 2^f assignments;
  (2) detailed balance pi(x) K(x,y) = pi(y) K(y,x) for every pair of states reachable;
  (3) a mutant (acceptance without the proposal terms) violates detailed balance.
"""
import itertools
from fractions import Fraction as F
from math import factorial

W = {0: F(1, 2), 1: F(3, 10)}          # weights of the two components A=0, B=1 (restricted to them)
ALPHA = F(1, 2)                         # Dirichlet-like rising-factorial marginal


def rising(a, n):
    r = F(1)
    for i in range(n):
        r *= a + i
    return r


def m(block_items, data):
    """Marginal of a block of categorical items: product over categories of rising factorials."""
    cnt = {}
    for i in block_items:
        cnt[data[i]] = cnt.get(data[i], 0) + 1
    n = len(block_items)
    r = F(1) / rising(3 * ALPHA, n)
    for c, k in cnt.items():
        r *= rising(ALPHA, k)
    return r


def target(assign, data, items):
    p = F(1)
    for c in (0, 1):
        blk = [i for i in items if assign[i] == c]
        p *= W[c] ** len(blk) * m(blk, data)
    return p


def proposal(path, order, data, anchors, free_order):
    """Sequential allocation density of `path` (dict item->comp) visiting free items in order."""
    blocks = {0: [anchors[0]], 1: [anchors[1]]}
    q = F(1)
    for it in order:
        sc = {}
        for c in (0, 1):
            sc[c] = W[c] * m(blocks[c] + [it], data) / m(blocks[c], data)
        z = sc[0] + sc[1]
        q *= sc[path[it]] / z
        blocks[path[it]].append(it)
    return q


def check(data, perm_average=True):
    n = len(data)
    items = list(range(n))
    anchors = (0, 1)
    free = items[2:]
    states = []
    for bits in itertools.product((0, 1), repeat=len(free)):
        a = {0: 0, 1: 1}
        a.update(dict(zip(free, bits)))
        states.append(a)
    orders = list(itertools.permutations(free))
    pi = [target(s, data, items) for s in states]
    Z = sum(pi)
    pi = [p / Z for p in pi]
    K = [[F(0)] * len(states) for _ in states]
    for xi, x in enumerate(states):
        for order in orders:
            pr_order = F(1, len(orders))
            for yi, y in enumerate(states):
                qxy = proposal(y, order, data, anchors, free)
                qyx = proposal(x, order, data, anchors, free)
                ratio = (pi[yi] * qyx) / (pi[xi] * qxy)
                acc = min(F(1), ratio)
                K[xi][yi] += pr_order * qxy * acc
        K[xi][xi] += 1 - sum(K[xi])
    for xi in range(len(states)):
        assert sum(K[xi]) == 1
        for yi in range(len(states)):
            assert pi[xi] * K[xi][yi] == pi[yi] * K[yi][xi], "detailed balance fails"
    for yi in range(len(states)):
        assert sum(pi[xi] * K[xi][yi] for xi in range(len(states))) == pi[yi]
    # (1) proposal normalised
    for order in orders:
        assert sum(proposal(y, order, data, anchors, free) for y in states) == 1
    # (4) mutant: drop proposal terms from the ratio
    bad = False
    Km = [[F(0)] * len(states) for _ in states]
    for xi in range(len(states)):
        for order in orders:
            for yi in range(len(states)):
                qxy = proposal(states[yi], order, data, anchors, free)
                Km[xi][yi] += F(1, len(orders)) * qxy * min(F(1), pi[yi] / pi[xi])
        Km[xi][xi] += 1 - sum(Km[xi])
    for xi in range(len(states)):
        for yi in range(len(states)):
            if pi[xi] * Km[xi][yi] != pi[yi] * Km[yi][xi]:
                bad = True
    assert bad, "mutant not detected"
    return len(states), len(orders)


if __name__ == "__main__":
    for data in ([0, 1, 0, 1, 2], [0, 0, 1, 2, 2, 1], [2, 2, 2, 0, 1, 0, 1]):
        ns, no = check(data)
        print(f"PASS data {data}: {ns} states, {no} visiting orders; detailed balance, invariance, "
              f"normalised proposal exact; mutant (no proposal terms) caught")
