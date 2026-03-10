"""
run_x12_x_1.sage — Compute A_f for y^2 = x^12 + x + 1 using pyrforest,
with N = 2^16 = 65536, and record the elapsed time.

Usage:
    sage pyrforest/run_x12_x_1.sage
"""

import timeit as timeit_module
timer = timeit_module.default_timer

import os as _os
load(_os.path.abspath("pyrforest/P7-avg-poly.sage"))

R.<x> = PolynomialRing(Integers())
f = x^12 + x + 1
C = HyperellipticCurve(f)
N = 2^16  # 65536

print(f"Computing A_f for y^2 = x^12 + x + 1 up to N = {N} using pyrforest ...")

start = timer()
answer = compute_A_f_avg_poly_from_curve(C, N)
elapsed = timer() - start

print(f"\nDone.  Elapsed time: {elapsed:.3f} seconds")
print(f"Number of primes processed: {len(answer)}")

# Print first 10 and last 10 results as a sanity check
keys = sorted(answer.keys())
print(f"\nFirst 10 primes:")
for p in keys[:10]:
    print(f"  p = {p}: {answer[p].list()}")
print(f"\nLast 10 primes:")
for p in keys[-10:]:
    print(f"  p = {p}: {answer[p].list()}")
