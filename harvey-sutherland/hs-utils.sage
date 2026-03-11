"""
hs-utils.sage — Utilities for the Harvey-Sutherland (arXiv:1410.5222) algorithm.

Usage:
    load("harvey-sutherland/hs-utils.sage")
"""

# -------------------------------------------------------------------------
# Parsing f into h
# -------------------------------------------------------------------------

def parse_f(f_coeffs):
    """
    Given coefficients of f as a list [f_0, f_1, ..., f_d] (index = degree),
    return (h_coeffs, r, c, e) where:
      h(x) = f(x) / x^c = h_0 + h_1*x + ... + h_r*x^r
      c = lowest nonzero degree of f (in {0, 1})
      r = d - c
      e = 2 - c (in {1, 2})
    """
    d = len(f_coeffs) - 1
    c = 0
    while c <= d and f_coeffs[c] == 0:
        c += 1
    r = d - c
    e = 2 - c
    h_coeffs = list(f_coeffs[c:])
    return h_coeffs, r, c, e


# -------------------------------------------------------------------------
# M^n_k polynomial matrix (Harvey-Sutherland eq. 3)
# -------------------------------------------------------------------------

def construct_M_k_poly(r, h_coeffs, P):
    """
    Construct the r x r polynomial matrix for M_k (Harvey-Sutherland eq. 7),
    where x = P.gen() stands in for k.

    M_k is the p-independent simplification of M^n_k, satisfying
    2*M^n_k == M_k (mod p).  It is this matrix -- not M^n_k -- that is
    fed into the remainder forest.

    M_k entries:
      - Subdiagonal (row i, col i-1): 2*k*h_0   --> 2*x*h_0
      - Last column (row i, col r-1): (r-i-2*k)*h_{r-i}
                                      --> (r-i)*h_{r-i} - 2*x*h_{r-i}
    """
    x = P.gen()
    M = matrix(P, r, r)

    for i in range(1, r):
        M[i, i-1] = 2 * x * h_coeffs[0]

    for i in range(r):
        M[i, r-1] = (r - i - 2*x) * h_coeffs[r - i]

    return M


def construct_M_k_n_poly(n, r, h_coeffs, P):
    """
    Construct the r x r polynomial matrix M^n_k (Harvey-Sutherland eq. 3),
    where x = P.gen() stands in for the loop variable k, and n is a fixed
    parameter (analogous to i in construct_T_bar_poly).

    M^n_k entries:
      - Subdiagonal (row i, col i-1): k * h_0          --> x * h_0
      - Last column (row i, col r-1): ((n+1)*(r-i)-k) * h_{r-i}
                                      --> ((n+1)*(r-i) - x) * h_{r-i}

    Evaluating at x=k recovers the concrete matrix M^n_k for that step.
    """
    x = P.gen()
    M = matrix(P, r, r)

    for i in range(1, r):
        M[i, i-1] = x * h_coeffs[0]

    for i in range(r):
        M[i, r-1] = ((n + 1) * (r - i) - x) * h_coeffs[r - i]

    return M
