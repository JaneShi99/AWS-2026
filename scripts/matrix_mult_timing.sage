"""
matrix_mult_timing.sage — Benchmark Sage integer matrix multiplication
for square matrices with large random integer entries.

Usage (run from repo root or anywhere):
    sage scripts/matrix_mult_timing.sage <dim_min> <dim_max> <n_trials> <bit_size>

Arguments:
    dim_min   Minimum matrix dimension (integer >= 1).
    dim_max   Maximum matrix dimension (integer >= dim_min).
    n_trials  Number of independent multiplications to time per dimension;
              results are averaged.  Use 1 for a quick run, 3–5 for more
              stable estimates.
    bit_size  Approximate bit-length of each matrix entry.  Entries are
              drawn uniformly at random from [0, 2^bit_size).

Notes:
    - Matrices are over ZZ (arbitrary-precision integers).
    - Timing uses time.process_time_ns() and is reported in milliseconds.
    - A power-law fit T(d) ≈ A·d^α is printed and overlaid on the plot.
    - The plot is saved to scripts/matrix_mult_timing.png.
    - For very large bit_size (e.g. 10^7) and dim >= 5, individual
      multiplications can take tens of seconds; plan accordingly.

Examples:
    sage scripts/matrix_mult_timing.sage 1 10 3 10000000
    sage scripts/matrix_mult_timing.sage 1 20 2 1000000
    sage scripts/matrix_mult_timing.sage 1 5  1 100000000
"""

import sys
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

USAGE = (
    "Usage: sage scripts/matrix_mult_timing.sage "
    "<dim_min> <dim_max> <n_trials> <bit_size>\n"
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
    DIM_MIN  = int(args[0])
    DIM_MAX  = int(args[1])
    N_TRIALS = int(args[2])
    BIT_SIZE = int(args[3])
except ValueError:
    print(f"Error: all arguments must be integers.\n{USAGE}", file=sys.stderr)
    sys.exit(1)

errors = []
if DIM_MIN < 1:
    errors.append("dim_min must be >= 1")
if DIM_MAX < DIM_MIN:
    errors.append("dim_max must be >= dim_min")
if N_TRIALS < 1:
    errors.append("n_trials must be >= 1")
if BIT_SIZE < 1:
    errors.append("bit_size must be >= 1")
if errors:
    for e in errors:
        print(f"Error: {e}", file=sys.stderr)
    print(USAGE, file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Benchmarking
# ---------------------------------------------------------------------------

def random_integer_matrix(d, bit_size):
    """Return a random d×d matrix over ZZ with entries in [0, 2^bit_size)."""
    bound = ZZ(2) ** bit_size
    return matrix(ZZ, d, d, [ZZ.random_element(bound) for _ in range(d * d)])


def time_mult(d, bit_size):
    """
    Time one multiplication of two fresh random d×d integer matrices.
    Returns elapsed process time in nanoseconds.
    """
    A = random_integer_matrix(d, bit_size)
    B = random_integer_matrix(d, bit_size)
    t0 = time.process_time_ns()
    _ = A * B
    t1 = time.process_time_ns()
    return t1 - t0


print(f"Benchmarking ZZ matrix multiplication")
print(f"  dimensions : d = {DIM_MIN} .. {DIM_MAX}")
print(f"  trials/dim : {N_TRIALS}")
print(f"  entry size : ~{BIT_SIZE:,} bits  (~{int(BIT_SIZE * 0.301):,} decimal digits)")
print()
sys.stdout.flush()

dims       = []
mean_times = []   # nanoseconds

for d in range(DIM_MIN, DIM_MAX + 1):
    times_ns = [time_mult(d, BIT_SIZE) for _ in range(N_TRIALS)]
    mean_ns  = sum(times_ns) / len(times_ns)
    dims.append(d)
    mean_times.append(mean_ns)
    trial_str = "  ".join(f"{t/1e6:.1f}" for t in times_ns)
    print(f"  d={d:3d}  mean={mean_ns/1e6:10.2f} ms    trials: [{trial_str}] ms")
    sys.stdout.flush()

# ---------------------------------------------------------------------------
# Power-law fit:  T(d) = A · d^α  via linear regression on log-log data
# ---------------------------------------------------------------------------

def fit_power_law(xs, ys):
    log_x = np.log(np.array(xs, dtype=float))
    log_y = np.log(np.array(ys, dtype=float))
    alpha, log_A = np.polyfit(log_x, log_y, 1)
    return np.exp(log_A), alpha


A, alpha = fit_power_law(dims, mean_times)
print(f"\nPower-law fit  T(d) ≈ A · d^α  (log-log linear regression):")
print(f"  α = {alpha:.3f}")
print(f"  (naive matrix multiplication is O(d^3); fast algorithms reach O(d^2.37))")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

d_arr    = np.array(dims, dtype=float)
times_ms = np.array(mean_times) / 1e6
d_fine   = np.linspace(DIM_MIN, DIM_MAX, 500)
fit_ms   = A * d_fine**alpha / 1e6

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle(
    rf"Sage $\mathbb{{Z}}$ matrix multiplication: process time vs dimension"
    f"\n(~{BIT_SIZE:,}-bit entries,  {N_TRIALS} trial(s)/dimension)",
    fontsize=12,
)

# Left: linear axes
ax = axes[0]
ax.scatter(dims, times_ms, color="firebrick", s=50, zorder=3, label="Measured")
ax.plot(d_fine, fit_ms, "--", color="firebrick", linewidth=1.5, alpha=0.8,
        label=rf"Fit: $T \propto d^{{{alpha:.2f}}}$")
ax.set_xlabel("Dimension $d$", fontsize=11)
ax.set_ylabel("Process time (ms)", fontsize=11)
ax.set_title("Linear scale", fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Right: log-log axes
ax = axes[1]
ax.scatter(d_arr, times_ms, color="firebrick", s=50, zorder=3, label="Measured")
ax.plot(d_fine, fit_ms, "--", color="firebrick", linewidth=1.5, alpha=0.8,
        label=rf"Fit: $T \propto d^{{{alpha:.2f}}}$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Dimension $d$  (log scale)", fontsize=11)
ax.set_ylabel("Process time (ms)  (log scale)", fontsize=11)
ax.set_title("Log-log scale", fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, which="both")

plt.tight_layout()

outfile = "scripts/matrix_mult_timing.png"
plt.savefig(outfile, dpi=150)
print(f"\nPlot saved to {outfile}")
