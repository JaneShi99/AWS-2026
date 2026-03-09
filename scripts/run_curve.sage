load("good-code/P7-utils.sage")
load("good-code/P7-avg-poly.sage")

R.<x> = ZZ[]
f = 3*x^4 + 6*x^3 + 6*x^2 + 5*x + 4
C = HyperellipticCurve(f)
print(f"Curve: y^2 = {f}")
print(f"Genus: {C.genus()}")
print()

# f.list() gives ascending coefficients: [const, x, x^2, ..., leading]
F_coeffs = list(f)
d = f.degree()

N = 100
result = compute_A_f_avg_poly(F_coeffs, d, N)

print(f"Hasse-Witt matrices A_f[p] for primes 7 <= p <= {N}:")
print()
for p in sorted(result.keys()):
    print(f"  p = {p}:  A_f = {result[p].list()}")
