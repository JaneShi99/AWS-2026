"""
profile_pyrforest.sage — cProfile a single compute_A_f_avg_poly_from_curve
call using the pyrforest implementation.

Usage (run from repo root):
    sage scripts/profile_pyrforest.sage <g> <N_max>
"""

import sys
import cProfile
import pstats
import io

args = sys.argv[1:]
if len(args) != 2 or args[0] in ("-h", "--help"):
    print("Usage: sage scripts/profile_pyrforest.sage <g> <N_max>")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

G    = int(args[0])
N    = int(args[1])

load("tests/hyperell_suite/random_generation/random_curves.sage")
load("pyrforest/P7-avg-poly.sage")

# Generate a good curve (retry if bad prime)
p_check = next_prime(max(G, 5))
print(f"Generating genus-{G} even-degree curve over GF({p_check})...")
while True:
    C = random_hyperelliptic_even_degree(p_check, G)
    result_check = compute_A_f_avg_poly_from_curve(C, p_check)
    if p_check in result_check:
        break

print(f"Profiling compute_A_f_avg_poly_from_curve(C, {N})...\n")

pr = cProfile.Profile()
pr.enable()
result = compute_A_f_avg_poly_from_curve(C, N)
pr.disable()

s = io.StringIO()
ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
ps.print_stats(30)
print(s.getvalue())
