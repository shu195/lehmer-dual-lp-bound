"""
GPU-Accelerated Self-Certifying Dual-LP Solver for Lehmer's Totient Problem (M >= 3)
====================================================================================
Uses Diagonal-Preconditioned PDHG (Chambolle-Pock / cuPDLP) in PyTorch (CUDA/float64)
to solve the Prime Divisibility DAG Dual LP at massive cutoffs (X = 2M, 5M, 10M+):

    min_{y >= 0, z >= 0}  sum_{e in E_X} y_e + sum_{p in V_X} z_p
    subject to            z_p + sum_{e ni p} y_e >= w_p   for all p in V_X

Every iterate y >= 0 is automatically an unconditional dual certificate once
z_p = max(0, w_p - (B y)_p). After GPU optimization, the best dual vector y_best
is verified on CPU in 50-digit directed-rounding Decimal + exact integer arithmetic
scaled by 10^12 (zero IEEE-754 float trust).

Usage (Colab / Kaggle / Local NVIDIA GPU):
    python lehmer_gpu_pdhg.py --X 5000000 --iters 15000 --M 3
"""
import argparse
import math
import time
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING

import numpy as np
import torch

SCALE = 10**12


def build_prime_dag_fast(X: int):
    """O(X log log X) SPF sieve to build V_X = {p prime : 5 <= p <= X} and E_X = {(q, p) : q | (p-1)}."""
    t0 = time.time()
    spf = np.arange(X + 1, dtype=np.int32)
    limit = int(math.isqrt(X))
    for i in range(2, limit + 1):
        if spf[i] == i:
            spf[i * i :: i] = np.minimum(spf[i * i :: i], i)

    is_prime = (spf == np.arange(X + 1, dtype=np.int32))
    is_prime[:5] = False
    ps = np.nonzero(is_prime)[0].astype(np.int64)

    prime_to_idx = np.full(X + 1, -1, dtype=np.int32)
    prime_to_idx[ps] = np.arange(len(ps), dtype=np.int32)

    u_list = []
    v_list = []
    for idx_p, p in enumerate(ps):
        m = int(p) - 1
        while m > 1:
            q = int(spf[m])
            if q >= 5:
                u = int(prime_to_idx[q])
                if u >= 0:
                    u_list.append(u)
                    v_list.append(idx_p)
            while m % q == 0:
                m //= q

    u_arr = np.array(u_list, dtype=np.int64)
    v_arr = np.array(v_list, dtype=np.int64)
    print(f"[DAG] X={X:,} | |V|={len(ps):,} primes | |E|={len(u_arr):,} edges | built in {time.time()-t0:.2f}s")
    return ps, u_arr, v_arr


@torch.no_grad()
def solve_dual_lp_pdhg_gpu(
    n_verts: int,
    u_np: np.ndarray,
    v_np: np.ndarray,
    w_np: np.ndarray,
    const_C: float,
    iters: int,
    eval_every: int,
    device: torch.device,
    label: str,
):
    """
    Edge-Scaled Diagonal-Preconditioned PDHG (cuPDLP style) + Water-Filling Init for:
        max_{x in [0, 1]^V} min_{y >= 0}  w^T x + y^T (1 - B^T x)
    with exact dual objective:
        F(y) = sum(y) + sum(max(0, w - B y))

    Key insight: w_p spans 9 orders of magnitude (10^-3 to 10^6), and deg(5) ~ pi(X)/4.
    By scaling each edge e = (u, v) (u < v) by s_e = max(w_v, 1e-3) and using exact
    Cauchy-Schwarz step sizes tau_i * sigma_e * (2 * row_sum_i) < 1, every dual variable
    y_hat_e = y_e / s_e lives in O(1) and converges in a few thousand GPU iterations!
    """
    t0 = time.time()
    u = torch.from_numpy(u_np).to(device=device, dtype=torch.long)
    v = torch.from_numpy(v_np).to(device=device, dtype=torch.long)
    w = torch.from_numpy(w_np).to(device=device, dtype=torch.float64)
    n_edges = u.numel()

    # 1. Exact closed-form GPU water-filling initialization:
    #    For every vertex v with in-degree d_in[v] > 0, distribute w[v] across its incoming edges.
    #    This immediately satisfies z_v = 0 for all v with incoming edges and credits small primes u!
    in_deg = torch.zeros(n_verts, device=device, dtype=torch.float64)
    ones_e = torch.ones(n_edges, device=device, dtype=torch.float64)
    in_deg.index_add_(0, v, ones_e)

    y = torch.clamp(w[v], min=0.0) / torch.clamp(in_deg[v], min=1.0)

    # 2. Edge-scaled diagonal preconditioning:
    #    Let s_e = max(w[v], 1e-2) be the natural scale of edge e = (u, v).
    #    Set sigma_e = eta * s_e, and tau_i = 0.98 / (eta * (sum_{e ni i} 2 * s_e)).
    #    Proof of convergence:
    #      <B y, x> = sum_e y_e (x_u + x_v)
    #               <= sum_e (|y_e| / sqrt(sigma_e)) * sqrt(sigma_e) * (|x_u| + |x_v|)
    #               <= ||y||_{sigma^{-1}} * sqrt( sum_e 2 * sigma_e * (x_u^2 + x_v^2) )
    #               =  ||y||_{sigma^{-1}} * sqrt( sum_i (sum_{e ni i} 2 * sigma_e) * x_i^2 )
    #               <= sqrt(0.98) * ||y||_{sigma^{-1}} * ||x||_{tau^{-1}}.
    s_e = torch.clamp(w[v], min=1e-2)
    row_sum = torch.zeros(n_verts, device=device, dtype=torch.float64)
    row_sum.index_add_(0, u, 2.0 * s_e)
    row_sum.index_add_(0, v, 2.0 * s_e)

    eta = 1.0
    sigma = eta * s_e
    tau = 0.98 / (eta * torch.clamp(row_sum, min=1e-12))

    # Initialize primal x from dual slack (vertices with positive slack want x_i = 1)
    By = torch.zeros(n_verts, device=device, dtype=torch.float64)
    By.index_add_(0, u, y)
    By.index_add_(0, v, y)
    x = (w > By).to(torch.float64)
    y_bar = y.clone()

    best_dual = float((y.sum() + torch.clamp(w - By, min=0.0).sum()).item())
    best_y = y.clone()
    print(f"  [{label}] init   0/{iters} | dual={best_dual:,.4f} | best_LB={const_C - best_dual:,.4f} | {time.time()-t0:.1f}s", flush=True)

    for it in range(1, iters + 1):
        # Primal update: x <- clip_[0,1](x + tau * (w - B y_bar))
        By.zero_()
        By.index_add_(0, u, y_bar)
        By.index_add_(0, v, y_bar)
        x.add_(tau * (w - By)).clamp_(0.0, 1.0)

        # Dual update: y_next = max(0, y + sigma * (x_u + x_v - 1))
        y_prev = y
        y = torch.clamp(y_prev + sigma * (x[u] + x[v] - 1.0), min=0.0)

        # Extrapolation: y_bar = 2 y - y_prev
        y_bar = 2.0 * y - y_prev

        if it % eval_every == 0 or it == iters:
            By.zero_()
            By.index_add_(0, u, y)
            By.index_add_(0, v, y)
            dual_obj = float((y.sum() + torch.clamp(w - By, min=0.0).sum()).item())
            if dual_obj < best_dual:
                best_dual = dual_obj
                best_y.copy_(y)
            lb_est = const_C - best_dual
            if it % (eval_every * 5) == 0 or it == iters:
                print(
                    f"  [{label}] iter {it:6d}/{iters} | dual={dual_obj:,.4f} | best_LB={lb_est:,.4f} | {time.time()-t0:.1f}s",
                    flush=True,
                )

    return best_y.cpu().numpy(), best_dual


def verify_exact_integer(ps, u_arr, v_arr, y_float, w_int_upper, const_C_int_lower, label: str):
    """50-digit directed-rounding + exact-integer verification (zero IEEE-754 float trust)."""
    y_int = [max(0, int(math.floor(float(val) * SCALE))) for val in y_float]
    aty = [0] * len(ps)
    for u, v, ye in zip(u_arr.tolist(), v_arr.tolist(), y_int):
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
    parser = argparse.ArgumentParser(description="GPU PDHG Dual-LP Certifier for Lehmer's Totient Problem")
    parser.add_argument("--X", type=int, default=2_000_000, help="Prime cutoff X (default: 2,000,000)")
    parser.add_argument("--M", type=int, default=3, help="Lower bound on M in n - 1 = M*phi(n) (default: 3)")
    parser.add_argument("--iters", type=int, default=15_000, help="PDHG iterations per objective (default: 15,000)")
    parser.add_argument("--eval-every", type=int, default=500, help="Evaluate exact dual objective every N iters")
    parser.add_argument("--skip-strict", action="store_true", help="Skip CPU 50-digit exact-integer verification")
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[DEVICE] Using {device}" + (f" ({torch.cuda.get_device_name(0)})" if device.type == "cuda" else ""))

    ps, u_arr, v_arr = build_prime_dag_fast(args.X)

    # Float64 weights for GPU PDHG
    c_f = math.log(args.X / (args.X - 1))
    lnM_f = math.log(args.M)
    R_f = math.log10(args.X) / c_f

    w_k_f = np.log(ps / (ps - 1.0)) / c_f - 1.0
    w_d_f = R_f * np.log(ps / (ps - 1.0)) - np.log10(ps)
    C_k_f = lnM_f / c_f
    C_d_f = R_f * lnM_f

    y_k, _ = solve_dual_lp_pdhg_gpu(len(ps), u_arr, v_arr, w_k_f, C_k_f, args.iters, args.eval_every, device, "omega(n)")
    y_d, _ = solve_dual_lp_pdhg_gpu(len(ps), u_arr, v_arr, w_d_f, C_d_f, args.iters, args.eval_every, device, "log10(n)")

    if args.skip_strict:
        return

    print("[VERIFY] Computing 50-digit directed-rounding exact-integer certificate on CPU...")
    with localcontext() as ctx:
        ctx.prec = 50
        ctx.rounding = ROUND_FLOOR
        lnM_lo = Decimal(args.M).ln()
        logX_lo = Decimal(args.X).log10()
        ctx.rounding = ROUND_CEILING
        c_hi = (Decimal(args.X) / Decimal(args.X - 1)).ln()

        ctx.rounding = ROUND_FLOOR
        C_k_int_lo = int(((lnM_lo / c_hi) * Decimal(SCALE)).to_integral_value(rounding=ROUND_FLOOR))
        R_lo = logX_lo / c_hi
        C_d_int_lo = int(((R_lo * lnM_lo) * Decimal(SCALE)).to_integral_value(rounding=ROUND_FLOOR))

        w_k_hi, w_d_hi = [], []
        for p in ps.tolist():
            ctx.rounding = ROUND_FLOOR
            logp_lo = Decimal(p).log10()
            ctx.rounding = ROUND_CEILING
            lnp_hi = (Decimal(p) / Decimal(p - 1)).ln()
            w_k_hi.append(int(((lnp_hi / c_hi - Decimal(1)) * Decimal(SCALE)).to_integral_value(rounding=ROUND_CEILING)))
            w_d_hi.append(int(((R_lo * lnp_hi - logp_lo) * Decimal(SCALE)).to_integral_value(rounding=ROUND_CEILING)))

    k_lb = verify_exact_integer(ps, u_arr, v_arr, y_k, w_k_hi, C_k_int_lo, f"Lehmer M>={args.M} omega(n)")
    d_lb = verify_exact_integer(ps, u_arr, v_arr, y_d, w_d_hi, C_d_int_lo, f"Lehmer M>={args.M} log10(n)")
    print(f"  => UNCONDITIONAL (M >= {args.M}): omega(n) >= {int(k_lb) + 1:,}, n > 10^{int(d_lb):,} ({int(d_lb) + 1:,} digits)")


if __name__ == "__main__":
    main()
