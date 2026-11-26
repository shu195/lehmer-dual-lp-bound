"""Generate 2 publication-quality figures (PDF + 300-DPI PNG) for the Lehmer M >= 3 amsart paper (X = 2.9M peak)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"C:\Users\Admin\.gemini\antigravity\scratch\lehmer_paper\figures"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 8.2,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
})

# Figure 1: Full Unimodal Saturation Curve across X in [0.321M, 7.0M]
X_vals = np.array([0.321, 0.6, 1.0, 2.0, 2.9, 3.0, 4.0, 5.0, 6.0, 7.0])
k_vals = np.array([26000, 39399, 53826, 77361, 85167, 85088, 79766, 67685, 50307, 28592])
d_vals = np.array([139.105, 220.630, 311.913, 466.953, 509.336, 508.786, 471.576, 385.914, 261.213, 103.831])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.6))

ax1.plot(X_vals, k_vals / 1000.0, marker="o", ms=4.5, lw=2.1, color="#1f4e79", label="Dual-LP Certified $\\omega(n)$")
ax1.scatter([2.9], [85.167], color="#d97706", s=55, zorder=5, label="Global Peak ($X=2.9$M)")
ax1.axhline(16.001, color="#b22222", ls="--", lw=1.2, label="Prior $M \\geq 3$ bound: 16,001")
ax1.annotate("85,167", xy=(2.9, 85.167), xytext=(4.1, 87.5),
             arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=1.0),
             fontweight="bold", color="#1f4e79")
ax1.set_title("(a) Certified Prime Factors $\\omega(n)$ vs. Cutoff $X$")
ax1.set_xlabel("DAG Cutoff $X$ (millions)")
ax1.set_ylabel("Certified Lower Bound (thousands)")
ax1.set_ylim(0, 98)
ax1.grid(alpha=0.25)
ax1.legend(loc="lower left", frameon=True)

ax2.plot(X_vals, d_vals, marker="s", ms=4.5, lw=2.1, color="#1b6b3a", label="Dual-LP Certified $\\log_{10}(n)$")
ax2.scatter([2.9], [509.336], color="#d97706", s=55, zorder=5, label="Global Peak ($X=2.9$M)")
ax2.axhline(8.171, color="#7f7f7f", ls=":", lw=1.2, label="Hagis (1988): 8,171 digits")
ax2.annotate("509,337 digits", xy=(2.9, 509.336), xytext=(4.0, 520),
             arrowprops=dict(arrowstyle="->", color="#1b6b3a", lw=1.0),
             fontweight="bold", color="#1b6b3a")
ax2.set_title("(b) Certified Decimal Digits $\\log_{10}(n)$ vs. Cutoff $X$")
ax2.set_xlabel("DAG Cutoff $X$ (millions)")
ax2.set_ylabel("Certified Digits (thousands)")
ax2.set_ylim(0, 580)
ax2.grid(alpha=0.25)
ax2.legend(loc="lower left", frameon=True)

fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig1_lehmer_convergence.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig1_lehmer_convergence.png"), bbox_inches="tight")
plt.close(fig)

# Figure 2: Comparison of M >= 3 Prime-Factor Records
labels = ["Hagis\n(1988)", "Prior Record\n($M \\geq 3$)", "Dual-LP\n($X=1$M)", "Dual-LP\n($X=2$M)", "Dual-LP Peak\n($X=2.9$M)"]
vals = [1991, 16001, 53826, 77361, 85167]
colors = ["#94a3b8", "#3b82f6", "#0ea5e9", "#10b981", "#15803d"]

fig, ax = plt.subplots(figsize=(7.6, 3.5))
bars = ax.bar(labels, vals, color=colors, width=0.55, edgecolor="#1e293b", linewidth=0.7)
for bar, v in zip(bars, vals):
    ax.text(bar.get_x() + bar.get_width()/2, v + 1800, f"{v:,}",
            ha="center", va="bottom", fontsize=9, fontweight="bold")
ax.set_ylabel("Minimum Distinct Prime Factors $\\omega(n)$")
ax.set_title("Lower Bounds on $\\omega(n)$ for Lehmer's Totient Problem ($M \\geq 3$)")
ax.set_ylim(0, 98000)
ax.grid(axis="y", alpha=0.25)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig2_lehmer_comparison.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig2_lehmer_comparison.png"), bbox_inches="tight")
plt.close(fig)
print("LEHMER 509K FIGURES GENERATED OK")
