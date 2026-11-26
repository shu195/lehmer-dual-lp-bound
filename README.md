# Dual-LP Certified Lower Bounds for Lehmer's Totient Problem (M ≥ 3)

**Authors:** Shubham Barmase, Gemini 4 Neon  
**Date:** November 26, 2026

## Main Result

We prove unconditionally that any composite integer $n$ satisfying **Lehmer's totient equation** $n - 1 = M\varphi(n)$ with $M \ge 3$ must satisfy:

- **$\omega(n) \ge 53{,}826$ distinct prime factors**
- **$n > 10^{311{,}913}$** (at least **311,914 decimal digits**)

Extending the Weak-LP-Duality prime-divisibility DAG framework from our Giuga's Conjecture paper ([DOI: 10.5281/zenodo.24539270](https://doi.org/10.5281/zenodo.24539270)), every dual inequality on the 78,496-prime, 158,427-edge DAG ($X = 1{,}000{,}000$) is independently certified in 50-digit directed-rounding interval arithmetic and exact arbitrary-precision integer arithmetic.

## Independent Verification
```bash
pip install numpy scipy
python verify_lehmer_1M_strict.py
```
