#!/usr/bin/env python3
import unittest
import pypaya2

class PypayaTestCase(unittest.TestCase):
    def assert_approx_equal(self, expected, actual):
        if abs(expected - actual) / abs(expected) > 1e-14:
            self.fail("expected {0} got {1}".format(expected, actual))

    def assert_sets_equal(self, expected, actual):
        expected = set(expected)
        actual = set(actual)
        extra = actual - expected
        missing = expected - actual
        if len(extra) > 0:
            self.fail("got extra elements: {}".format(extra))
        if len(missing) > 0:
            self.fail("missing the elements: {}".format(missing))
