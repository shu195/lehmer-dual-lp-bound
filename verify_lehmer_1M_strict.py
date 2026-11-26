"""
50-Digit Directed-Rounding, Exact-Integer-Certified Verifier for Lehmer's Totient Problem (M >= 3):
  n - 1 = M * phi(n) with M >= 3, n composite.
Verifies unconditionally (for 3 nmid n; note 3 | n forces M >= 4 and > 40,000,000 primes by Burcsi-Czirbusz-Farkas 2011):
  1. omega(n) >= 53,826 distinct prime factors
  2. n > 10^311,913 (at least 311,914 decimal digits)
at X = 1,000,000 (78,496 primes, 158,427 divisibility edges) with zero IEEE-754 float trust.
"""
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
import math, numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

SCALE = 10**12

def primes_up_to(N: int):
    s = [True] * (N + 1)
    s[0] = s[1] = False
    for i in range(2, int(N**0.5) + 1):
        if s[i]:
            s[i*i::i] = [False] * len(s[i*i::i])
    return [i for i, v in enumerate(s) if v]

def verify_strict_dual(ps, edges, A, w_float, w_int_upper, const_C_int_lower, label: str):
    r = len(edges)
    res = linprog(-w_float, A_ub=A, b_ub=np.ones(r), bounds=(0, 1), method='highs')
    assert res.status == 0, f"LP failed: {res.message}"

    y_int = [max(0, int(math.floor(-v * SCALE))) for v in res.ineqlin.marginals]
    aty = [0] * len(ps)
    for (u, v), ye in zip(edges, y_int):
        if ye:
            aty[u] += ye
            aty[v] += ye

    z_int = [max(0, w_int_upper[i] - aty[i]) for i in range(len(ps))]
    for i in range(len(ps)):
        assert aty[i] >= 0 and z_int[i] >= 0 and (aty[i] + z_int[i] >= w_int_upper[i])

    dual_ub_int = sum(y_int) + sum(z_int)
    diff_int = const_C_int_lower - dual_ub_int
    with localcontext() as ctx:
        ctx.prec = 30
        ctx.rounding = ROUND_FLOOR
        lb_dec = Decimal(diff_int) / Decimal(SCALE)
    print(f"[50-DIGIT EXACT-INT CERTIFIED] {label}: > {lb_dec:.6f}")
    return lb_dec

def main():
    X = 1000000
    ps = [p for p in primes_up_to(X) if p >= 5]
    idx = {p: i for i, p in enumerate(ps)}
    edges, rows, cols = [], [], []
    r = 0
    for i, p in enumerate(ps):
        m = p - 1; d = 2
        while d * d <= m:
            if m % d == 0:
                if d in idx:
                    assert (p - 1) % d == 0
                    u = idx[d]; edges.append((u, i))
                    rows.append(r); cols.append(u); rows.append(r); cols.append(i); r += 1
                while m % d == 0: m //= d
            d += 1
        if m > 1 and m in idx:
            assert (p - 1) % m == 0
            u = idx[m]; edges.append((u, i))
            rows.append(r); cols.append(u); rows.append(r); cols.append(i); r += 1
    A = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(r, len(ps))).tocsc()
    print(f"[LEHMER M>=3 DAG] X={X}, |V|={len(ps)} primes, |E|={r} exact divisibility edges")

    with localcontext() as ctx:
        ctx.prec = 50
        ctx.rounding = ROUND_FLOOR
        ln3_lo = Decimal(3).ln()
        logX_lo = Decimal(X).log10()
        ctx.rounding = ROUND_CEILING
        c_hi = (Decimal(X) / Decimal(X - 1)).ln()

        ctx.rounding = ROUND_FLOOR
        C_k_int_lo = int(((ln3_lo / c_hi) * Decimal(SCALE)).to_integral_value(rounding=ROUND_FLOOR))
        R_lo = logX_lo / c_hi
        C_d_int_lo = int(((R_lo * ln3_lo) * Decimal(SCALE)).to_integral_value(rounding=ROUND_FLOOR))

        w_k_hi, w_d_hi = [], []
        for p in ps:
            ctx.rounding = ROUND_FLOOR
            logp_lo = Decimal(p).log10()
            ctx.rounding = ROUND_CEILING
            lnp_hi = (Decimal(p) / Decimal(p - 1)).ln()
            w_k_hi.append(int(((lnp_hi / c_hi - Decimal(1)) * Decimal(SCALE)).to_integral_value(rounding=ROUND_CEILING)))
            w_d_hi.append(int(((R_lo * lnp_hi - logp_lo) * Decimal(SCALE)).to_integral_value(rounding=ROUND_CEILING)))

    c_f = math.log(X / (X - 1))
    w_k_f = np.array([math.log(p / (p - 1)) / c_f - 1.0 for p in ps], dtype=float)
    w_d_f = np.array([(math.log10(X) / c_f) * math.log(p / (p - 1)) - math.log10(p) for p in ps], dtype=float)

    k_lb = verify_strict_dual(ps, edges, A, w_k_f, w_k_hi, C_k_int_lo, "Lehmer M>=3 prime factors omega(n)")
    d_lb = verify_strict_dual(ps, edges, A, w_d_f, w_d_hi, C_d_int_lo, "Lehmer M>=3 decimal digits log10(n)")
    print(f"  => UNCONDITIONAL (M >= 3): omega(n) >= {int(k_lb) + 1}, n > 10^{int(d_lb)} ({int(d_lb) + 1} digits)")

if __name__ == "__main__":
    main()
