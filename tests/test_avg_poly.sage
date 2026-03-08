"""
test_avg_poly.sage — Tests for compute_A_f_avg_poly against hyperell_suite cases.

Run via: sage tests/runtests.sage  (from repo root)
"""

import unittest

load("tests/hyperell_suite/sage/loader.sage")
load("good-code/P7-utils.sage")
load("good-code/P7-avg-poly.sage")

CASE_FILE = "tests/hyperell_suite/cases/random_p17_to_50.json"


class TestAvgPolyRandomP17To50(unittest.TestCase):

    def test_cases(self):
        cases = load_cases(CASE_FILE)
        for case in cases:
            with self.subTest(curve=str(case.curve)):
                p = case.curve.base_ring().characteristic()
                result = compute_A_f_avg_poly_from_curve(case.curve, p)

                # p may be a bad prime (e.g. p divides f(0)), in which case
                # the algorithm skips it and result[p] won't exist
                if p not in result:
                    self.skipTest(f"p={p} is a bad prime for this curve")

                A_f = result[p]
                g = A_f.nrows()

                # Compute det(1 - T*A_f) over GF(p)
                R_poly.<T> = PolynomialRing(GF(p))
                M = identity_matrix(R_poly, g) - T * A_f.change_ring(R_poly)
                det_coeffs = M.det().list()
                computed = det_coeffs + [GF(p)(0)] * (g + 1 - len(det_coeffs))

                # L(T) ≡ det(1 - T*A_f) mod p, so compare the first g+1 coefficients
                expected = [c % p for c in case.expected_Lpoly[:g+1]]

                self.assertEqual(computed, expected)
