"""
complexity_vs_genus_pyrforest.sage — Benchmark the pyrforest implementation
of compute_A_f_avg_poly_from_curve across a range of genera, using random
even-degree (2g+2) hyperelliptic curves.

For each genus g we use the smallest prime p > max(g, 5) as p_check (to
guarantee p_check > g and p_check > 5).  A_f is computed for all primes up
to N_max.  Two random curves are generated per genus by default; if p_check
is a bad prime for a curve a new one is drawn automatically.

Measures process time in nanoseconds, plots time vs genus, and fits a
power-law T(g) ≈ A·g^α to estimate the complexity exponent.

Usage (run from repo root):
    sage scripts/complexity_vs_genus_pyrforest.sage <g_min> <g_max> <n_per_genus> <N_max>

Arguments:
    g_min        Minimum genus (integer >= 1).
    g_max        Maximum genus (integer >= g_min).
    n_per_genus  Number of random curves to benchmark per genus; timings are
                 averaged.  Use 1 for a quick run, 3–5 for more stable results.
    N_max        Upper bound passed to compute_A_f_avg_poly_from_curve as N.
                 The algorithm computes A_f for every prime p <= N_max.
                 Must be >= next_prime(max(g_max, 5)) so that at least one
                 valid prime (p > g) is computed for the largest genus.

Notes:
    - Curves are generated over GF(p) where p = next_prime(max(g, 5)), which
      guarantees p > g (the algorithm requires p > g) and p > 5.
    - If that p is a bad prime for a curve (p | f(0)), a new curve is drawn.
    - The plot is saved to scripts/complexity_vs_genus_pyrforest_<timestamp>.png.
    - A power-law fit T(g) ≈ A·g^α is printed and overlaid on the plot.

Examples:
    sage scripts/complexity_vs_genus_pyrforest.sage 1 30 2 41
    sage scripts/complexity_vs_genus_pyrforest.sage 1 100 2 101
    sage scripts/complexity_vs_genus_pyrforest.sage 10 50 3 200
"""

import sys
import time
import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

USAGE = (
    "Usage: sage scripts/complexity_vs_genus_pyrforest.sage "
    "<g_min> <g_max> <n_per_genus> <N_max>\n"
    "Run with --help for full documentation."
)

args = sys.argv[1:]

if args and args[0] in ("-h", "--help"):
    print(__doc__)
    sys.exit(0)

if len(args) != 4:
    print(f"Error: expected 4 arguments, got {len(args)}.\n{USAGE}", file=sys.stderr)
    sys.exit(1)

try:
    G_MIN       = int(args[0])
    G_MAX       = int(args[1])
    N_PER_GENUS = int(args[2])
    N_MAX       = int(args[3])
except ValueError:
    print(f"Error: all arguments must be integers.\n{USAGE}", file=sys.stderr)
    sys.exit(1)

errors = []
if G_MIN < 1:
    errors.append("g_min must be >= 1")
if G_MAX < G_MIN:
    errors.append("g_max must be >= g_min")
if N_PER_GENUS < 1:
    errors.append("n_per_genus must be >= 1")
if N_MAX < 1:
    errors.append("N_max must be >= 1")
else:
    min_valid_p = next_prime(max(G_MAX, 5))
    if N_MAX < min_valid_p:
        errors.append(
            f"N_max ({N_MAX}) must be >= next_prime(max(g_max, 5)) = {min_valid_p} "
            f"so that at least one valid prime p > g_max is included"
        )
if errors:
    for e in errors:
        print(f"Error: {e}", file=sys.stderr)
    print(USAGE, file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Load pyrforest implementation and curve generator
# ---------------------------------------------------------------------------

load("tests/hyperell_suite/random_generation/random_curves.sage")
load("pyrforest/P7-avg-poly.sage")

# ---------------------------------------------------------------------------
# Benchmarking
# ---------------------------------------------------------------------------

def benchmark_genus(g, N):
    """
    Generate N_PER_GENUS random even-degree curves of genus g, compute A_f
    for all primes up to N using the pyrforest implementation, and return a
    list of process-time measurements in nanoseconds.

    Curves are generated over GF(p_check) where p_check = next_prime(max(g, 5))
    to guarantee p_check > g (required by the algorithm).  If p_check turns
    out to be a bad prime for a given curve, a new curve is drawn.
    """
    p_check = next_prime(max(g, 5))
    times_ns = []

    for _ in range(N_PER_GENUS):
        while True:
            C = random_hyperelliptic_even_degree(p_check, g)
            t0 = time.process_time_ns()
            result = compute_A_f_avg_poly_from_curve(C, N)
            t1 = time.process_time_ns()
            if p_check in result:
                times_ns.append(t1 - t0)
                break

    return times_ns


print(f"Benchmarking pyrforest implementation, even-degree (2g+2) curves")
print(f"  genus range : g = {G_MIN} .. {G_MAX}")
print(f"  curves/genus: {N_PER_GENUS}")
print(f"  N_max       : {N_MAX}  (compute A_f for all primes up to N_max)")
print()
sys.stdout.flush()

genera     = []
mean_times = []   # nanoseconds

for g in range(G_MIN, G_MAX + 1):
    p_check = next_prime(max(g, 5))
    times_ns = benchmark_genus(g, N_MAX)
    mean_ns  = sum(times_ns) / len(times_ns)
    genera.append(g)
    mean_times.append(mean_ns)
    print(f"  g={g:3d}  p_check={p_check:4d}  mean={mean_ns/1e6:9.2f} ms  ({N_PER_GENUS} curves)")
    sys.stdout.flush()

# ---------------------------------------------------------------------------
# Power-law fit:  T(g) = A · g^α  via linear regression on log-log data
# ---------------------------------------------------------------------------

def fit_power_law(genera, times_ns):
    log_g = np.log(np.array(genera, dtype=float))
    log_t = np.log(np.array(times_ns, dtype=float))
    alpha, log_A = np.polyfit(log_g, log_t, 1)
    return np.exp(log_A), alpha


A, alpha = fit_power_law(genera, mean_times)
print(f"\nPower-law fit  T(g) ≈ A · g^α  (log-log linear regression):")
print(f"  α = {alpha:.3f}")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

g_arr    = np.array(genera, dtype=float)
times_ms = np.array(mean_times) / 1e6
g_fine   = np.linspace(G_MIN, G_MAX, 500)
fit_ms   = A * g_fine**alpha / 1e6

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle(
    rf"pyrforest $A_f$ algorithm: process time vs genus"
    f"\n(even degree $2g+2$,  $N={N_MAX}$,  {N_PER_GENUS} curve(s)/genus)",
    fontsize=12,
)

# Left: linear axes
ax = axes[0]
ax.scatter(genera, times_ms, color="darkorange", s=30, zorder=3, label="Measured")
ax.plot(g_fine, fit_ms, "--", color="darkorange", linewidth=1.5, alpha=0.8,
        label=rf"Fit: $T \propto g^{{{alpha:.2f}}}$")
ax.set_xlabel("Genus $g$", fontsize=11)
ax.set_ylabel("Process time (ms)", fontsize=11)
ax.set_title("Linear scale", fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Right: log-log axes
ax = axes[1]
ax.scatter(g_arr, times_ms, color="darkorange", s=30, zorder=3, label="Measured")
ax.plot(g_fine, fit_ms, "--", color="darkorange", linewidth=1.5, alpha=0.8,
        label=rf"Fit: $T \propto g^{{{alpha:.2f}}}$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Genus $g$  (log scale)", fontsize=11)
ax.set_ylabel("Process time (ms)  (log scale)", fontsize=11)
ax.set_title("Log-log scale", fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, which="both")

plt.tight_layout()

timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
outfile = f"scripts/complexity_vs_genus_pyrforest_{timestamp}.png"
plt.savefig(outfile, dpi=150)
print(f"\nPlot saved to {outfile}")
