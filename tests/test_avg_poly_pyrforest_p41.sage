"""
test_avg_poly_pyrforest_p41.sage — Tests for the pyrforest implementation of
compute_A_f_avg_poly_from_curve against high-genus cases at p=41 (genus 1
through 30, even degree only).

Run via: sage tests/runtests.sage  (from repo root)
"""

import unittest

load("tests/hyperell_suite/sage/loader.sage")
load("pyrforest/P7-avg-poly.sage")

EVEN_CASE_FILE = "tests/hyperell_suite/cases/hyperelliptic_p41_even_degree.json"


class TestAvgPolyPyrforestP41EvenDegree(unittest.TestCase):

    def test_cases(self):
        cases = load_cases(EVEN_CASE_FILE)
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

                R_poly.<T> = PolynomialRing(GF(p))
                M = identity_matrix(R_poly, g) - T * A_f.change_ring(R_poly)
                det_coeffs = M.det().list()
                computed = det_coeffs + [GF(p)(0)] * (g + 1 - len(det_coeffs))

                expected = [c % p for c in case.expected_Lpoly[:g+1]]

                self.assertEqual(computed, expected)
