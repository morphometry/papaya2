#!/usr/bin/env python3
import unittest
import numpy as np
import pypaya2

class PypayaTestCase(unittest.TestCase):
    def assert_approx_equal(self, expected, actual, accuracy=1e-14):
        expected = np.asarray(expected)
        actual = np.asarray(actual)
        delta = expected - actual
        if (abs(delta) > abs(expected * accuracy)).any():
            self.fail("expected {0} got {1} [delta = {2}]".format(expected, actual, delta))

    def assert_equal(self, expected, actual):
        self.assertEqual(expected, actual)

    def assert_sets_equal(self, expected, actual):
        expected = set(expected)
        actual = set(actual)
        extra = actual - expected
        missing = expected - actual
        if len(extra) > 0:
            self.fail("got extra elements: {}".format(extra))
        if len(missing) > 0:
            self.fail("missing the elements: {}".format(missing))
