"""
hs-matrix.sage — Algorithms ComputeHasseWittMatrix and ComputeHasseWittMatrices
from Harvey-Sutherland (arXiv:1410.5222, Section 6).

Usage (from repo root):
    sage harvey-sutherland/hs-matrix.sage
"""

import os as _os
if _os.path.exists("harvey-sutherland/hs-first-row.sage"):
    load(_os.path.abspath("harvey-sutherland/hs-first-row.sage"))
elif _os.path.exists("hs-first-row.sage"):
    load(_os.path.abspath("hs-first-row.sage"))


def _deduce_wp_from_first_rows(first_rows, shifts, p, g):
    """
    Lemma WpFromFirstRows: reconstruct W_p from first rows of g translated curves.

    first_rows[i] = W_p^1(shifts[i]) as a list of g elements of GF(p).
    shifts[i] are distinct mod p.

    For each column jp = 0,...,g-1 (paper column j = jp+1):
      (a) Compute w_{jp+1}(a_i) from the first jp columns already known
          using eq. (wj) from the paper (empty sum for jp=0).
      (b) Solve V(shifts) * col = [W_p^1(a_i)[jp] - w_{jp+1}(a_i)]_i.

    Returns W_p as a g x g matrix over GF(p).
    """
    Fp = GF(p)
    V = matrix(Fp, g, g, lambda i, k: Fp(shifts[i])^k)
    V_inv = V.inverse()
    W = matrix(Fp, g, g)

    for jp in range(g):  # paper column index j = jp+1
        # eq. (wj): w_{jp+1}(a) =
        #   sum_{k'=0}^{g-1} sum_{l'=0}^{jp-1} C(jp,l') * (-1)^{jp-l'} * a^{k'+jp-l'} * W[k',l']
        # (empty sum for jp=0)
        wj_vals = []
        for i in range(g):
            a = Fp(shifts[i])
            val = Fp(0)
            for kp in range(g):
                for lp in range(jp):
                    val += (Fp(binomial(jp, lp) * (-1)^(jp - lp))
                            * a^(kp + jp - lp)
                            * W[kp, lp])
            wj_vals.append(val)

        rhs = vector(Fp, [first_rows[i][jp] - wj_vals[i] for i in range(g)])
        W[:, jp] = V_inv * rhs

    return W


def compute_hasse_witt_matrix(f_bar_coeffs, p):
    """
    Algorithm ComputeHasseWittMatrix (Harvey-Sutherland §6, before Thm 6.1).

    Given integer coefficients of f_bar = f mod p and a prime p >= g,
    computes the full Hasse-Witt matrix W_p of y^2 = f_bar(x) over GF(p).

    Steps:
      1. Choose shifts a_i = 0,...,g-1 (distinct mod p since p >= g).
      2. Compute W_p^1(a_i) via ComputeHasseWittFirstRow on y^2 = f_bar(x+a_i).
      3. Reconstruct W_p via Lemma WpFromFirstRows.

    Returns W_p as a g x g matrix over GF(p).
    """
    Fp = GF(p)
    Rx = PolynomialRing(Fp, 'x')
    x = Rx.gen()

    f_bar_fp = Rx([Fp(c) for c in f_bar_coeffs])
    d = f_bar_fp.degree()
    g = (d - 1) // 2

    assert p >= g, f"p={p} must be >= g={g}"

    shifts = list(range(g))

    first_rows = []
    for a in shifts:
        f_shifted = f_bar_fp(x + Fp(a))
        f_shifted_coeffs = [Integer(c) for c in f_shifted.list()]
        first_rows.append(compute_hasse_witt_first_row(f_shifted_coeffs, p))

    return _deduce_wp_from_first_rows(first_rows, shifts, p, g)


def compute_hasse_witt_matrices(f_coeffs, N):
    """
    Algorithm ComputeHasseWittMatrices (Harvey-Sutherland §6).

    Given integer coefficients of f and N, computes W_p for all good primes p <= N.

    Steps:
      1. For odd primes p < g of good reduction, compute W_p directly:
           W_p[i,j] = [x^{p(i+1)-(j+1)}] f^{(p-1)/2} mod p  (0-indexed i,j).
      2. Choose shifts a_i = 0,...,g-1.  Let S be the set of primes p in [g,N]
         of good reduction that divide some nonzero integer f(a_i).
         (Differences a_i - a_j are < g <= p, so they contribute nothing to S.)
      3. For p in S, compute W_p via ComputeHasseWittMatrix (single-prime).
      4. For p in [g,N] of good reduction not in S, run ComputeHasseWittFirstRows
         once per shift on f(x+a_i) over ZZ, then deduce W_p via Lemma WpFromFirstRows.
         If a shift is inadmissible for some p (unexpected), fall back to step 3.

    Returns {p: W_p as a g x g matrix over GF(p)} for all good primes p <= N.
    """
    d = len(f_coeffs) - 1
    g = (d - 1) // 2

    # Good reduction: f mod p is squarefree (computed on demand)
    def is_good(p):
        if p == 2:
            return False
        Fp = GF(p)
        Rx = PolynomialRing(Fp, 'x')
        f_p = Rx([Fp(c) for c in f_coeffs])
        return gcd(f_p, f_p.derivative()).degree() == 0

    result = {}

    # ------------------------------------------------------------------
    # Step 1: small primes 3 <= p < g (direct Hasse-Witt definition)
    # ------------------------------------------------------------------
    for p in prime_range(3, g):
        if not is_good(p):
            continue
        Fp = GF(p)
        Rx = PolynomialRing(Fp, 'x')
        f_p = Rx([Fp(c) for c in f_coeffs])
        f_power = f_p^((p - 1) // 2)
        W = matrix(Fp, g, g)
        for i in range(g):
            for j in range(g):
                idx = p * (i + 1) - (j + 1)
                if idx >= 0:
                    W[i, j] = f_power[idx]
        result[p] = W

    # Shifts a_i = 0,...,g-1 (differences < g <= p, so never contribute to S)
    shifts = list(range(g))

    # ------------------------------------------------------------------
    # Step 2: identify S = primes in [g,N] dividing some nonzero f(a_i)
    # ------------------------------------------------------------------
    R_ZZ = PolynomialRing(ZZ, 't')
    f_ZZ = R_ZZ(f_coeffs)
    f_at_shifts = [ZZ(f_ZZ(a)) for a in shifts]

    s_product = ZZ(1)
    for val in f_at_shifts:
        if val != 0:
            s_product = lcm(s_product, abs(val))

    S = {p for p in prime_range(g, N + 1) if s_product % p == 0}

    # ------------------------------------------------------------------
    # Step 3: primes in S — single-prime algorithm
    # ------------------------------------------------------------------
    for p in sorted(S):
        if not is_good(p):
            continue
        result[p] = compute_hasse_witt_matrix(f_coeffs, p)

    # ------------------------------------------------------------------
    # Step 4: remaining primes in [g,N] — forest algorithm + Vandermonde
    # ------------------------------------------------------------------
    t = R_ZZ.gen()

    # Run ComputeHasseWittFirstRows once per shift for all primes up to N
    first_rows_by_shift = []
    for a in shifts:
        f_shifted_coeffs = [Integer(c) for c in f_ZZ(t + a).list()]
        forest = compute_hasse_witt_first_rows(f_shifted_coeffs, N)
        first_rows_by_shift.append(forest)

    for p in prime_range(g, N + 1):
        if p in S or p in result:
            continue
        if not is_good(p):
            continue

        # For p not in S: each f(a_i) is coprime to p, so p is admissible
        # for every shifted curve and must appear in every forest result.
        if not all(p in first_rows_by_shift[i] for i in range(g)):
            # Unexpected: fall back to single-prime algorithm
            result[p] = compute_hasse_witt_matrix(f_coeffs, p)
            continue

        first_rows = [first_rows_by_shift[i][p] for i in range(g)]
        result[p] = _deduce_wp_from_first_rows(first_rows, shifts, p, g)

    return result


# ----------------------------------------------------------------------
# Demo
# ----------------------------------------------------------------------
import os as _os
if 'hs-matrix' in _os.path.basename(sys.argv[0]):
    R.<t> = PolynomialRing(ZZ)
    f = -(t^12 - t^10 + 6*t^9 - 7*t^8 + 5*t^7 + t^6 - t^5 + t^4 - t^3 + t^2 - t + 1)
    f_coeffs = [Integer(c) for c in f.list()]

    print("=== ComputeHasseWittMatrix (single prime) ===")
    for p in [11, 13, 17, 97]:
        W = compute_hasse_witt_matrix(f_coeffs, p)
        print(f"p={p}:")
        print(W)
        print()

    print("=== ComputeHasseWittMatrices (all primes up to N=30) ===")
    result = compute_hasse_witt_matrices(f_coeffs, 30)
    for p, W in sorted(result.items()):
        print(f"p={p}:")
        print(W)
        print()
