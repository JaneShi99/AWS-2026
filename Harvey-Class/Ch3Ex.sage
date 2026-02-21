p = 101

def count_points_naive(q, f):
    """
    Count points on y^2 = f(x) over F_q using naive enumeration.
    
    For each x in F_q:
    - If f(x) = 0: contributes 1 point (x, 0)
    - If f(x) is a nonzero square: contributes 2 points (x, ±sqrt(f(x)))
    - If f(x) is not a square: contributes 0 points
    
    Plus points at infinity: for degree 7 (odd), there's 1 point at infinity.
    """
    Fq = GF(q)
    
    # Helper to count points for each x value
    def pts(x):
        fx = f(Fq(x))
        return 1 if fx == 0 else (2 if fx.is_square() else 0)
    
    # Affine points using list comprehension + point at infinity
    return sum(pts(x) for x in Fq) + 1

# Define the polynomial ring and the curve
R.<x> = PolynomialRing(ZZ)

# From y^2 + x^7 - x^6 + 6x^5 - 7x^4 + 5x^3 + x^2 - x + 1 = 0
# We get y^2 = -(x^7 - x^6 + 6x^5 - 7x^4 + 5x^3 + x^2 - x + 1)
f = -(x^7 - x^6 + 6*x^5 - 7*x^4 + 5*x^3 + x^2 - x + 1)

print("="*60)
print("Point Counting on Genus 3 Hyperelliptic Curve")
print("="*60)
print(f"\nCurve: y^2 = {f}")
print(f"Base prime: p = {p}")
print(f"Genus: g = 3")
print()

# Count points over F_p, F_p^2, F_p^3
for k in range(1, 4):
    q = p^k
    N = count_points_naive(q, f)
    print(f"#C(F_{{{p}}}^{k}) = #C(F_{{{q}}}) = {N}")

print()
print("="*60)
print("Verification using SageMath's built-in HyperellipticCurve")
print("="*60)

# Alternative: Use SageMath's built-in hyperelliptic curve functionality
for k in range(1, 4):
    q = p^k
    Fq = GF(q)
    Rq.<X> = PolynomialRing(Fq)
    
    # For HyperellipticCurve(f, h), we have y^2 + h*y = f
    # Our equation is y^2 = -(x^7 - x^6 + 6x^5 - 7x^4 + 5x^3 + x^2 - x + 1)
    # So h = 0 and f = -(X^7 - X^6 + 6*X^5 - 7*X^4 + 5*X^3 + X^2 - X + 1)
    fq = -(X^7 - X^6 + 6*X^5 - 7*X^4 + 5*X^3 + X^2 - X + 1)
    
    try:
        C = HyperellipticCurve(fq)
        N_sage = C.count_points(1)[0]  # count_points returns list for extensions
        print(f"#C(F_{{{q}}}) = {N_sage} (SageMath HyperellipticCurve)")
    except Exception as e:
        print(f"F_{{{q}}}: Error using HyperellipticCurve: {e}")
