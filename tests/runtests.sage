"""
runtests.sage — Top-level test runner.

Usage (from repo root):
    sage tests/runtests.sage
"""

import unittest

load("tests/test_avg_poly.sage")

suite = unittest.TestLoader().loadTestsFromTestCase(TestAvgPolyRandomP17To50)
unittest.TextTestRunner(verbosity=2).run(suite)
