# Timing Experiment: Point Counting on Hyperelliptic Curves
# Measuring runtime as a function of p (prime), g (genus), r (extension degree)

import time

def random_hyperelliptic_poly(Fq, g):
    """
    Generate a random monic polynomial of degree 2g+1 for a hyperelliptic curve y^2 = f(x).
    Ensures f has no repeated roots (curve is nonsingular).
    """
    R = PolynomialRing(Fq, 'x')
    x = R.gen()
    while True:
        # Random monic polynomial of degree 2g+1
        f = x^(2*g + 1) + R.random_element(degree=2*g)
        # Check it's squarefree (no repeated roots)
        if f.gcd(f.derivative()) == 1:
            return f

def count_points_sage(f):
    """
    Count points on y^2 = f(x) using Sage's builtin HyperellipticCurve.count_points().
    Uses more efficient algorithms than naive enumeration.
    """
    C = HyperellipticCurve(f)
    return C.count_points(1)[0]

def time_point_count(p, g, r, num_trials=3):
    """
    Time the point counting for a random hyperelliptic curve.
    Returns average time over num_trials.
    """
    q = p^r
    Fq = GF(q)
    
    times = []
    for _ in range(num_trials):
        f = random_hyperelliptic_poly(Fq, g)
        
        start = time.time()
        N = count_points_sage(f)
        elapsed = time.time() - start
        times.append(elapsed)
    
    return sum(times) / len(times)

# ============================================================
# EXPERIMENT 1: Vary p (prime), fix g and r
# ============================================================
def experiment_vary_p(primes, g_fixed, r_fixed):
    """Vary p while keeping g and r fixed."""
    print(f"\n{'='*60}")
    print(f"EXPERIMENT 1: Varying p (prime)")
    print(f"Fixed: g = {g_fixed}, r = {r_fixed}")
    print(f"{'='*60}")
    print(f"{'p':>10} | {'q = p^r':>15} | {'Time (s)':>12} | {'#C(F_q)':>15}")
    print("-" * 60)
    
    results = []
    for p in primes:
        q = p^r_fixed
        Fq = GF(q)
        f = random_hyperelliptic_poly(Fq, g_fixed)
        
        start = time.time()
        N = count_points_sage(f)
        elapsed = time.time() - start
        
        print(f"{p:>10} | {q:>15} | {elapsed:>12.6f} | {N:>15}")
        results.append((p, q, elapsed, N))
    
    return results

# ============================================================
# EXPERIMENT 2: Vary g (genus), fix p and r
# ============================================================
def experiment_vary_g(p_fixed, genera, r_fixed):
    """Vary g while keeping p and r fixed."""
    print(f"\n{'='*60}")
    print(f"EXPERIMENT 2: Varying g (genus)")
    print(f"Fixed: p = {p_fixed}, r = {r_fixed}")
    print(f"{'='*60}")
    print(f"{'g':>10} | {'deg(f)=2g+1':>12} | {'Time (s)':>12} | {'#C(F_q)':>15}")
    print("-" * 60)
    
    q = p_fixed^r_fixed
    Fq = GF(q)
    
    results = []
    for g in genera:
        f = random_hyperelliptic_poly(Fq, g)
        
        start = time.time()
        N = count_points_sage(f)
        elapsed = time.time() - start
        
        print(f"{g:>10} | {2*g+1:>12} | {elapsed:>12.6f} | {N:>15}")
        results.append((g, 2*g+1, elapsed, N))
    
    return results

# ============================================================
# EXPERIMENT 3: Vary r (extension degree), fix p and g
# ============================================================
def experiment_vary_r(p_fixed, g_fixed, extensions):
    """Vary r while keeping p and g fixed."""
    print(f"\n{'='*60}")
    print(f"EXPERIMENT 3: Varying r (extension degree)")
    print(f"Fixed: p = {p_fixed}, g = {g_fixed}")
    print(f"{'='*60}")
    print(f"{'r':>10} | {'q = p^r':>15} | {'Time (s)':>12} | {'#C(F_q)':>15}")
    print("-" * 60)
    
    results = []
    for r in extensions:
        q = p_fixed^r
        Fq = GF(q)
        f = random_hyperelliptic_poly(Fq, g_fixed)
        
        start = time.time()
        N = count_points_sage(f)
        elapsed = time.time() - start
        
        print(f"{r:>10} | {q:>15} | {elapsed:>12.6f} | {N:>15}")
        results.append((r, q, elapsed, N))
    
    return results

# ============================================================
# RUN ALL EXPERIMENTS
# ============================================================

print("="*60)
print("HYPERELLIPTIC CURVE POINT COUNTING: TIMING EXPERIMENTS")
print("Algorithm: Sage's HyperellipticCurve.count_points()")
print("="*60)

# Experiment 1: Vary p from small to larger primes
# Fix g=2 (genus 2), r=1 (base field)
primes_to_test = [101, 251, 503, 1009, 2003, 4001, 8009]
results_p = experiment_vary_p(primes_to_test, g_fixed=2, r_fixed=1)

# Experiment 2: Vary g from 1 to 10
# Fix p=101, r=1
genera_to_test = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
results_g = experiment_vary_g(p_fixed=101, genera=genera_to_test, r_fixed=1)

# Experiment 3: Vary r from 1 to 4
# Fix p=101, g=2
# Note: r=4 means q = 101^4 ≈ 100 million, which will be slow!
extensions_to_test = [1, 2, 3]  # Keeping small to avoid very long runtimes
results_r = experiment_vary_r(p_fixed=101, g_fixed=2, extensions=extensions_to_test)

# ============================================================
# ANALYSIS: Expected complexity
# ============================================================
print(f"\n{'='*60}")
print("COMPLEXITY ANALYSIS")
print("="*60)
print("""
Sage uses sophisticated algorithms for point counting:
  - For small fields: enumeration
  - For larger fields: Schoof-Elkies-Atkin style algorithms (for elliptic curves)
  - For hyperelliptic curves: various methods including p-adic algorithms

Expected scaling with Sage's algorithms:
  - Varying p: May scale better than O(p) for large p
  - Varying g: Complexity increases with genus (more complex Jacobian)
  - Varying r: Extension fields are more expensive
""")

# ============================================================
# SUMMARY TABLE
# ============================================================
print(f"\n{'='*60}")
print("SUMMARY: Timing Results")
print("="*60)

print("\nExperiment 1 (vary p, g=2, r=1):")
for p, q, t, N in results_p:
    print(f"  p={p:>5}: {t:.4f}s")

print("\nExperiment 2 (vary g, p=101, r=1):")
for g, deg, t, N in results_g:
    print(f"  g={g:>2} (deg={deg:>2}): {t:.4f}s")

print("\nExperiment 3 (vary r, p=101, g=2):")
for r, q, t, N in results_r:
    print(f"  r={r}: q={q:>10}, {t:.4f}s")
