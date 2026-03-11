"""
hs-first-row.sage — Algorithm ComputeHasseWittFirstRows from Harvey-Sutherland
(arXiv:1410.5222, Section 4).

Usage (from repo root):
    sage harvey-sutherland/hs-first-row.sage
"""

from pyrforest import remainder_forest

import os as _os
if _os.path.exists("harvey-sutherland/hs-utils.sage"):
    load(_os.path.abspath("harvey-sutherland/hs-utils.sage"))
elif _os.path.exists("hs-utils.sage"):
    load(_os.path.abspath("hs-utils.sage"))


def compute_hasse_witt_first_rows(f_coeffs, N):
    """
    Algorithm ComputeHasseWittFirstRows (Harvey-Sutherland §4).

    Given integer coefficients of f and a bound N, computes W_p^1 (the first
    row of the Hasse-Witt matrix) for all admissible primes p <= N
    simultaneously using the accumulating remainder forest.

    A prime p is admissible if it is odd, does not divide h_0 (the constant
    term of h = f/x^c), and is a prime of good reduction.

    Returns a dict {p: [w_{11}, ..., w_{1g}]} with entries in GF(p).

    Steps follow the paper directly:
      1. Set up moduli m_n = p if p=2n+1 admissible, else 1.
      2. Run remainder forest to get u_n = V_0 M'_1 ... M'_n mod p.
      3. Run remainder forest (or use Wilson) to get delta_n = (en)! mod p.
      4. Scale u_n and extract last g entries in reverse order.
    """
    h_coeffs, r, c, e = parse_f(f_coeffs)
    d = len(f_coeffs) - 1
    g = (d - 1) // 2
    h0 = h_coeffs[0]

    P.<x> = ZZ[]

    # Polynomial matrix for M_k (eq. 7) -- independent of p
    M_poly = construct_M_k_poly(r, h_coeffs, P)

    # V_0 = [0, ..., 0, 1]  as a 1 x r matrix
    V0 = matrix(ZZ, 1, r)
    V0[0, r - 1] = 1

    # Admissible primes: odd primes p <= N not dividing h_0
    admissible = [p for p in prime_range(3, N + 1) if gcd(ZZ(h0), p) == 1]

    # ------------------------------------------------------------------
    # Step 2: u_n = V_0 * M'_1 * ... * M'_n  mod p
    #
    # e=1: M'_k = M_k, product k = 1 ... n = (p-1)/2
    #      --> remainder_forest range(1, (p+1)//2)
    # e=2: M'_k = M_{2k-1} * M_{2k}; equivalently use M_k individually
    #      for k = 1 ... 2n = p-1
    #      --> remainder_forest range(1, p)
    # ------------------------------------------------------------------
    if e == 1:
        k_upper_u = lambda p: (p + 1) // 2
    else:
        k_upper_u = lambda p: p

    forest_u = remainder_forest(
        M_poly,
        lambda p: p,
        k_upper_u,
        kbase=1,
        indices=admissible,
        V=V0,
    )

    # ------------------------------------------------------------------
    # Step 3: delta_n = (en)! mod p
    #
    # e=2: (p-1)! = -1 mod p  by Wilson's theorem
    # e=1: n! mod p via a scalar remainder forest with matrix [[x]]
    #      (product of 1*2*...*n = n!)
    # ------------------------------------------------------------------
    if e == 2:
        delta = {p: -1 for p in admissible}
    else:
        M_fact = matrix(P, 1, 1, [x])
        V_fact = matrix(ZZ, 1, 1, [1])
        forest_delta = remainder_forest(
            M_fact,
            lambda p: p,
            lambda p: (p + 1) // 2,
            kbase=1,
            indices=admissible,
            V=V_fact,
        )
        delta = {p: Integer(forest_delta[p][0, 0]) for p in admissible}

    # ------------------------------------------------------------------
    # Step 4: for each admissible prime p = 2n+1, output the last g
    # entries of  (2|p)^e / ((h0|p)^(e-1) * delta_n) * u_n  mod p
    # in reverse order.
    # ------------------------------------------------------------------
    result = {}
    for p in admissible:
        d_n = delta[p] % p
        if d_n == 0:
            continue

        leg_2  = kronecker(2,       p)
        leg_h0 = kronecker(ZZ(h0),  p)

        # scalar = (2|p)^e / ((h0|p)^(e-1) * delta_n)  mod p
        numerator   = Integer(leg_2**e) % p
        denominator = (Integer(leg_h0**(e - 1)) * d_n) % p
        if denominator == 0:
            continue
        scalar = (numerator * inverse_mod(int(denominator), p)) % p

        # Scale u_n (a 1 x r matrix) and extract last g entries reversed
        u_n = forest_u[p]
        Fp  = GF(p)
        scaled = [Fp(scalar) * Fp(Integer(u_n[0, j])) for j in range(r)]
        W_p1 = list(reversed(scaled[r - g:]))
        result[p] = W_p1

    return result


def compute_hasse_witt_first_row(f_bar_coeffs, p):
    """
    Algorithm ComputeHasseWittFirstRow (Harvey-Sutherland §4, before Thm 4.5).

    Naively computes W_p^1 for a single prime p by iterating the recurrence
    directly over k = 1 ... e_bar*n.  Runs in O(g*p) time.

    f_bar_coeffs: integer coefficients of f_bar = f mod p (as a list,
                  index = degree), representing the curve y^2 = f_bar(x)
                  over F_p.  The leading and constant nonzero coefficients
                  must be nonzero mod p.
    p:            the prime.

    Returns W_p^1 as a list of g elements of GF(p).
    """
    Fp = GF(p)
    f_bar_fp = [Fp(c) for c in f_bar_coeffs]

    # Derive h_bar, r_bar, c_bar, e_bar from f_bar over F_p
    d_bar = len(f_bar_fp) - 1
    c_bar = 0
    while c_bar <= d_bar and f_bar_fp[c_bar] == 0:
        c_bar += 1
    r_bar = d_bar - c_bar
    e_bar = 2 - c_bar
    h_bar = f_bar_fp[c_bar:]   # h_bar[i] = coefficient of x^i in h_bar(x)

    g   = (d_bar - 1) // 2
    n   = (p - 1) // 2
    h0  = h_bar[0]

    # Step 1: u_0 = [0, ..., 0, 1],  delta_0 = 1
    u     = [Fp(0)] * r_bar
    u[-1] = Fp(1)
    delta = Fp(1)

    # Step 2: for k = 1 ... e_bar*n  do  u_k = u_{k-1} * M_bar_k,  delta_k *= k
    P_fp = PolynomialRing(Fp, 'x_fp')
    M_k_poly_fp = construct_M_k_poly(r_bar, h_bar, P_fp)
    u_mat = matrix(Fp, 1, r_bar, u)
    for k in range(1, e_bar * n + 1):
        k_fp = Fp(k)
        M_k = M_k_poly_fp.apply_map(lambda poly: poly(k_fp))
        u_mat = u_mat * M_k
        delta = delta * k_fp
    u = list(u_mat[0])

    # Sparse custom matmul (exploits structure of M_bar_k):
    # M_bar_k is sparse: subdiagonal entries = 2k*h0,
    #                    last-column entry at row i = (r-i-2k)*h_{r-i}.
    # So  (u * M_bar_k)[j] = u[j+1] * 2k*h0          for j < r-1
    #     (u * M_bar_k)[r-1] = sum_i  u[i] * (r-i-2k) * h_{r-i}
    # for k in range(1, e_bar * n + 1):
    #     k_fp      = Fp(k)
    #     two_k_h0  = Fp(2) * k_fp * h0
    #
    #     new_u = [Fp(0)] * r_bar
    #
    #     for j in range(r_bar - 1):
    #         new_u[j] = u[j + 1] * two_k_h0
    #
    #     last = Fp(0)
    #     for i in range(r_bar):
    #         last += u[i] * Fp(r_bar - i - 2*k) * h_bar[r_bar - i]
    #     new_u[r_bar - 1] = last
    #
    #     u     = new_u
    #     delta = delta * k_fp

    # Step 3: scale and extract last g entries in reverse order
    leg_2  = kronecker(2,                    p)
    leg_h0 = kronecker(Integer(f_bar_fp[c_bar]), p)

    scalar = Fp(leg_2**e_bar) / (Fp(leg_h0**(e_bar - 1)) * delta)

    scaled = [scalar * v for v in u]
    return list(reversed(scaled[r_bar - g:]))


# ----------------------------------------------------------------------
# Demo
# ----------------------------------------------------------------------
import os as _os
if 'hs-first-row' in _os.path.basename(sys.argv[0]):
    R.<t> = PolynomialRing(ZZ)
    f = -(t^12 - t^10 + 6*t^9 - 7*t^8 + 5*t^7 + t^6 - t^5 + t^4 - t^3 + t^2 - t + 1)
    f_coeffs = [Integer(c) for c in f.list()]

    N = 100
    result = compute_hasse_witt_first_rows(f_coeffs, N)
    for p, row in sorted(result.items()):
        print(f"p={p}: W_p^1 = {row}")
