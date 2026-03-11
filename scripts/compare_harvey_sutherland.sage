"""
compare_harvey_sutherland.sage — Benchmark and compare the pyrforest
(remainder forest) implementation of compute_A_f_avg_poly and the
Harvey-Sutherland implementation of compute_hasse_witt_matrices across
a range of genera, using random even-degree (2g+2) hyperelliptic curves.

The same random curves are used for both implementations at each genus so
the comparison is fair.  A power-law T(g) ≈ A·g^α is fitted to each series
and both are overlaid on a single plot.

Usage (run from repo root):
    sage scripts/compare_harvey_sutherland.sage <g_min> <g_max> <n_per_genus> <N_max>

Arguments:
    g_min        Minimum genus (integer >= 1).
    g_max        Maximum genus (integer >= g_min).
    n_per_genus  Number of random curves to benchmark per genus; timings are
                 averaged.  Use 1 for a quick run, 3–5 for more stable results.
    N_max        Upper bound passed to both algorithms as N.  Both compute the
                 Hasse-Witt matrix for every prime p <= N_max.
                 Must be >= next_prime(max(g_max, 5)) so that at least one
                 valid prime (p > g) is included for the largest genus.

Notes:
    - Curves are generated over GF(p) where p = next_prime(max(g, 5)),
      guaranteeing p > g and p > 5.  If p is bad for a curve a new one is drawn.
    - The same curve instances are timed against both implementations.
    - The plot is saved to scripts/compare_harvey_sutherland_<timestamp>.png.

Examples:
    sage scripts/compare_harvey_sutherland.sage 1 10 2 41
    sage scripts/compare_harvey_sutherland.sage 1 20 2 101
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
    "Usage: sage scripts/compare_harvey_sutherland.sage "
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
# Load both implementations
# ---------------------------------------------------------------------------

load("tests/hyperell_suite/random_generation/random_curves.sage")
load("pyrforest/P7-avg-poly.sage")
pyrforest_compute = compute_A_f_avg_poly_from_curve

load("harvey-sutherland/hs-matrix.sage")

def hs_compute(C, N):
    f_poly, _ = C.hyperelliptic_polynomials()
    f_coeffs = [Integer(c) for c in f_poly]
    return compute_hasse_witt_matrices(f_coeffs, N)

# ---------------------------------------------------------------------------
# Benchmarking
# ---------------------------------------------------------------------------

def generate_good_curves(g):
    """
    Generate N_PER_GENUS random even-degree curves of genus g that are valid
    for both implementations at p_check = next_prime(max(g, 5)).
    Returns a list of HyperellipticCurve objects.
    """
    p_check = next_prime(max(g, 5))
    curves = []
    while len(curves) < N_PER_GENUS:
        C = random_hyperelliptic_even_degree(p_check, g)
        if p_check in pyrforest_compute(C, p_check) and p_check in hs_compute(C, p_check):
            curves.append(C)
    return curves


def time_impl(impl_fn, curves, N):
    """
    Time impl_fn(C, N) on each curve.  Returns list of process-time measurements
    in nanoseconds (only for runs where p_check appears in the result).
    """
    p_check = next_prime(max(curves[0].genus(), 5))
    times_ns = []
    for C in curves:
        t0 = time.process_time_ns()
        result = impl_fn(C, N)
        t1 = time.process_time_ns()
        if p_check in result:
            times_ns.append(t1 - t0)
    return times_ns


print(f"Comparing pyrforest vs Harvey-Sutherland implementations")
print(f"  genus range : g = {G_MIN} .. {G_MAX}")
print(f"  curves/genus: {N_PER_GENUS}  (same curves used for both)")
print(f"  N_max       : {N_MAX}")
print()
print(f"  {'g':>3}  {'p_check':>7}  {'pyrforest (ms)':>14}  {'harvey-suth (ms)':>16}")
print(f"  {'-'*3}  {'-'*7}  {'-'*14}  {'-'*16}")
sys.stdout.flush()

genera        = []
pf_mean_times = []   # nanoseconds
hs_mean_times = []   # nanoseconds

for g in range(G_MIN, G_MAX + 1):
    p_check = next_prime(max(g, 5))
    curves  = generate_good_curves(g)

    pf_times = time_impl(pyrforest_compute, curves, N_MAX)
    hs_times = time_impl(hs_compute,        curves, N_MAX)

    pf_mean = sum(pf_times) / len(pf_times)
    hs_mean = sum(hs_times) / len(hs_times)

    genera.append(g)
    pf_mean_times.append(pf_mean)
    hs_mean_times.append(hs_mean)

    print(f"  g={g:3d}  p={p_check:4d}  {pf_mean/1e6:14.2f}  {hs_mean/1e6:16.2f}")
    sys.stdout.flush()

# ---------------------------------------------------------------------------
# Power-law fit:  T(g) = A · g^α  via linear regression on log-log data
# ---------------------------------------------------------------------------

def fit_power_law(genera, times_ns):
    log_g = np.log(np.array(genera, dtype=float))
    log_t = np.log(np.array(times_ns, dtype=float))
    alpha, log_A = np.polyfit(log_g, log_t, 1)
    return np.exp(log_A), alpha


pf_A, pf_alpha = fit_power_law(genera, pf_mean_times)
hs_A, hs_alpha = fit_power_law(genera, hs_mean_times)

print(f"\nPower-law fit  T(g) ≈ A · g^α  (log-log linear regression):")
print(f"  pyrforest      : α = {pf_alpha:.3f}")
print(f"  harvey-suth    : α = {hs_alpha:.3f}")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

g_arr   = np.array(genera, dtype=float)
pf_ms   = np.array(pf_mean_times) / 1e6
hs_ms   = np.array(hs_mean_times) / 1e6
g_fine  = np.linspace(G_MIN, G_MAX, 500)
pf_fit  = pf_A * g_fine**pf_alpha / 1e6
hs_fit  = hs_A * g_fine**hs_alpha / 1e6

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle(
    rf"Hasse-Witt: pyrforest $A_f$ vs Harvey-Sutherland"
    f"\n(even degree $2g+2$,  $N={N_MAX}$,  {N_PER_GENUS} curve(s)/genus)",
    fontsize=12,
)

for ax, log_scale in [(axes[0], False), (axes[1], True)]:
    ax.scatter(g_arr, pf_ms, color="darkorange", s=30, zorder=3,
               label="pyrforest")
    ax.plot(g_fine, pf_fit, "--", color="darkorange", linewidth=1.5, alpha=0.8,
            label=rf"pyrforest fit: $\alpha={pf_alpha:.2f}$")

    ax.scatter(g_arr, hs_ms, color="steelblue", s=30, zorder=3,
               label="harvey-sutherland")
    ax.plot(g_fine, hs_fit, "--", color="steelblue", linewidth=1.5, alpha=0.8,
            label=rf"harvey-suth fit: $\alpha={hs_alpha:.2f}$")

    ax.set_xlabel("Genus $g$" + ("  (log scale)" if log_scale else ""), fontsize=11)
    ax.set_ylabel("Process time (ms)" + ("  (log scale)" if log_scale else ""), fontsize=11)
    ax.set_title("Log-log scale" if log_scale else "Linear scale", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which="both" if log_scale else "major")
    if log_scale:
        ax.set_xscale("log")
        ax.set_yscale("log")

plt.tight_layout()

timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
outfile = f"scripts/compare_harvey_sutherland_{timestamp}.png"
plt.savefig(outfile, dpi=150)
print(f"\nPlot saved to {outfile}")
