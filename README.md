# Dual-LP Certified Lower Bounds for Lehmer's Totient Problem (M ≥ 3)

**Authors:** Shubham Barmase, Gemini 4 Neon  
**Date:** November 26, 2026

## Main Result

We prove unconditionally that any composite integer $n$ satisfying **Lehmer's totient equation** $n - 1 = M\varphi(n)$ with $M \ge 3$ must satisfy:

- **$\omega(n) \ge 85{,}167$ distinct prime factors**
- **$n > 10^{509{,}336}$** (at least **509,337 decimal digits**)

Extending the Weak-LP-Duality prime-divisibility DAG framework from our Giuga's Conjecture paper ([DOI: 10.5281/zenodo.24539270](https://doi.org/10.5281/zenodo.24539270)), we map the complete saturation curve up to $X = 10{,}000{,}000$ via GPU-accelerated edge-scaled PDHG and certify every dual inequality at the global peak $X = 2{,}900{,}000$ ($210{,}107$ primes, $442{,}070$ edges) in 50-digit directed-rounding interval arithmetic and exact arbitrary-precision integer arithmetic.

## Independent Verification
```bash
# GPU-accelerated PDHG + 50-digit exact-integer verification at X = 2,900,000 (~8s on Tesla T4)
pip install numpy torch
python lehmer_gpu_pdhg.py --X 2900000 --iters 5000 --eval-every 1000 --M 3

# CPU HiGHS + 50-digit exact-integer verification at X = 1,000,000
pip install numpy scipy
python verify_lehmer_1M_strict.py
```
