"""
p_rank.sage — Wrappers for computing the p-rank of a hyperelliptic curve.

The p-rank of a curve C / F_p is computed as:

    p_rank(C, p) = rank( A_f(p)^g )

where A_f(p) is the g×g Hasse-Witt matrix of C at p, and the rank is taken
over F_p (i.e. as a matrix over GF(p)).

IMPORTANT: This formula is only valid for the base field F_p. It does NOT
generalise to extension fields F_{p^r} with r > 1; in that setting the
p-rank must be read off from the Newton polygon of the characteristic
polynomial of Frobenius.

Two functions are provided:
  - p_rank_naive(C, p)       computes the p-rank for a single prime p using the
                              O(p^{1/2})-per-prime sqrt-time algorithm
                              (compute_A_f_fast from P6-Af.sage).

  - p_rank_avg_poly(C, N)    uses the average-polynomial-time algorithm from
                              P7-avg-poly.sage (compute_A_f_avg_poly_from_curve)
                              and returns a dict {p: p_rank} for all good primes
                              up to N.
"""

import os as _os

# Load the sqrt-time single-prime implementation.
if _os.path.exists("P6-Af.sage"):
    load(_os.path.abspath("P6-Af.sage"))
else:
    load(_os.path.abspath("good-code/P6-Af.sage"))

# Load the average-polynomial-time multi-prime implementation.
# (This also loads P7-utils.sage internally.)
if _os.path.exists("P7-avg-poly.sage"):
    load(_os.path.abspath("P7-avg-poly.sage"))
else:
    load(_os.path.abspath("good-code/P7-avg-poly.sage"))


def _p_rank_from_A_f(A_f_mod_p, p):
    """
    Given the Hasse-Witt matrix A_f already reduced mod p (as a Sage matrix),
    return the p-rank: rank(A_f^g) over GF(p).

    NOTE: The formula rank(A_f^g) is only valid over the base field F_p,
    not over extension fields F_{p^r}.
    """
    g = A_f_mod_p.nrows()
    Fp = GF(p)
    A = A_f_mod_p.change_ring(Fp)
    return (A ** g).rank()


def p_rank_naive(C, p):
    """
    Compute the p-rank of the hyperelliptic curve C at the prime p using the
    O(p^{1/2})-per-prime sqrt-time algorithm (compute_A_f_fast from P6-Af.sage).

    Parameters
    ----------
    C : HyperellipticCurve
        A hyperelliptic curve defined over a finite field or ZZ.
    p : int
        A prime > 5 that is not a bad prime for C (i.e. p does not divide the
        constant term of the defining polynomial f).

    Returns
    -------
    int
        The p-rank of C, an integer in {0, 1, ..., g}.
    """
    F_coeffs_poly, _ = C.hyperelliptic_polynomials()
    F_coeffs = [ZZ(c) for c in F_coeffs_poly]
    A_f = compute_A_f_fast(F_coeffs, p)
    return _p_rank_from_A_f(A_f, p)


def p_rank_avg_poly(C, N):
    """
    Compute the p-rank of the hyperelliptic curve C for every good prime p <= N
    using the average-polynomial-time algorithm (compute_A_f_avg_poly_from_curve
    from P7-avg-poly.sage).

    Parameters
    ----------
    C : HyperellipticCurve
        A hyperelliptic curve defined over a finite field or ZZ.
    N : int
        Upper bound: p-ranks are computed for all good primes p <= N.

    Returns
    -------
    dict
        A dictionary {p: p_rank} for each good prime p <= N for which the
        algorithm produced a Hasse-Witt matrix.
    """
    A_f_dict = compute_A_f_avg_poly_from_curve(C, N)
    return {p: _p_rank_from_A_f(A_f, p) for p, A_f in A_f_dict.items()}
