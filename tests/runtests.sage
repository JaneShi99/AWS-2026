"""
runtests.sage — Top-level test runner.

Usage (from repo root):
    sage tests/runtests.sage
"""

import unittest

load("tests/test_avg_poly.sage")
load("tests/test_avg_poly_p41.sage")
load("tests/test_avg_poly_pyrforest_p41.sage")

suite = unittest.TestSuite([
    unittest.TestLoader().loadTestsFromTestCase(TestAvgPolyRandomP17To50),
    unittest.TestLoader().loadTestsFromTestCase(TestAvgPolyP41OddDegree),
    unittest.TestLoader().loadTestsFromTestCase(TestAvgPolyP41EvenDegree),
    unittest.TestLoader().loadTestsFromTestCase(TestAvgPolyPyrforestP41EvenDegree),
])
unittest.TextTestRunner(verbosity=2).run(suite)
