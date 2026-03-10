"""
fit_genus_data.sage — Fit hardcoded benchmark data to four models:
  1. T = A * g^2          (fixed power, OLS)
  2. T = A * g^3          (fixed power, OLS)
  3. T = A * g^alpha      (free power law, log-log linear regression)
  4. T = A * exp(b * g)   (exponential, log-y linear regression)

Data from: complexity_vs_genus_pyrforest, g=2..5, N_max=2000000, 1 curve/genus.

Usage (run from repo root):
    sage scripts/fit_genus_data.sage
"""

import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

genera    = np.array([2, 3, 4, 5], dtype=float)
times_ms  = np.array([1926250.45, 2497712.90, 3412248.75, 5142391.10])

# ---------------------------------------------------------------------------
# Fits: T(g) = A * g^k, least-squares A = sum(T*g^k) / sum(g^(2k))
# ---------------------------------------------------------------------------

def fit_fixed_power(genera, times, k):
    """Fit T = A * g^k via least squares, returning A."""
    gk = genera ** k
    A = np.dot(times, gk) / np.dot(gk, gk)
    return A

A2 = fit_fixed_power(genera, times_ms, 2)
A3 = fit_fixed_power(genera, times_ms, 3)

# Free power law: log(T) = log(A) + alpha*log(g)
log_g = np.log(genera)
log_t = np.log(times_ms)
alpha, log_A_pow = np.polyfit(log_g, log_t, 1)
A_pow = np.exp(log_A_pow)

# Exponential: log(T) = log(A) + b*g
b, log_A_exp = np.polyfit(genera, log_t, 1)
A_exp = np.exp(log_A_exp)

print(f"Fit  T(g) = A * g^2:         A = {A2:.2f} ms")
print(f"Fit  T(g) = A * g^3:         A = {A3:.2f} ms")
print(f"Fit  T(g) = A * g^alpha:     A = {A_pow:.2f} ms,  alpha = {alpha:.3f}")
print(f"Fit  T(g) = A * exp(b*g):    A = {A_exp:.2f} ms,  b = {b:.3f}")

# Residual sum of squares (in original ms, not log space)
def rss(pred):
    r = times_ms - pred
    return np.dot(r, r)

print(f"\nResidual sum of squares (ms^2):")
print(f"  g^2 fit:         {rss(A2 * genera**2):.4e}")
print(f"  g^3 fit:         {rss(A3 * genera**3):.4e}")
print(f"  free power law:  {rss(A_pow * genera**alpha):.4e}")
print(f"  exponential:     {rss(A_exp * np.exp(b * genera)):.4e}")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

g_fine      = np.linspace(1.8, 5.5, 300)
fit2_ms     = A2 * g_fine**2
fit3_ms     = A3 * g_fine**3
fit_pow_ms  = A_pow * g_fine**alpha
fit_exp_ms  = A_exp * np.exp(b * g_fine)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle(
    r"pyrforest $A_f$: process time vs genus  ($N=2{,}000{,}000$, 1 curve/genus)",
    fontsize=12,
)

for ax, yscale in zip(axes, ["linear", "log"]):
    ax.scatter(genera, times_ms, color="darkorange", s=60, zorder=5, label="Measured")
    ax.plot(g_fine, fit2_ms, "--",  color="steelblue",  linewidth=1.8,
            label=rf"$T = {A2:.0f}\,g^2$")
    ax.plot(g_fine, fit3_ms, "-.",  color="seagreen",   linewidth=1.8,
            label=rf"$T = {A3:.0f}\,g^3$")
    ax.plot(g_fine, fit_pow_ms, ":", color="crimson",   linewidth=2.0,
            label=rf"$T = {A_pow:.0f}\,g^{{{alpha:.2f}}}$  (free)")
    ax.plot(g_fine, fit_exp_ms, "-", color="purple",    linewidth=1.5, alpha=0.8,
            label=rf"$T = {A_exp:.0f}\,e^{{{b:.3f}\,g}}$")
    ax.set_xlabel("Genus $g$", fontsize=11)
    ax.set_ylabel("Process time (ms)", fontsize=11)
    ax.set_yscale(yscale)
    ax.set_title(f"{'Linear' if yscale == 'linear' else 'Log-y'} scale", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

plt.tight_layout()

timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
outfile = f"scripts/fit_genus_data_{timestamp}.png"
plt.savefig(outfile, dpi=150)
print(f"\nPlot saved to {outfile}")
