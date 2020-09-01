#!/usr/bin/env python3
import unittest
import numpy as np
import pypaya2

class PolygonTest(unittest.TestCase):
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

    vertices = [
        [2., 3.], [.25, .5], [-3., 4.], [-2., 0.], [2., 0.]
    ]

    expected_area = 67. / 8

    expected_imt = [
        18.9509878231867863325879476083,
        0.,
        4.03586622195516544754949115540 + 3.83613097404496685735871370308j,
        -2.25242470620683687651763840592 + 7.91410261263648605530287871331j,
        2.23638454508277526919840618339 + 2.16586100869600377036856029938j,
        0.70071569649771790984433739720 - 10.30885074737037023789259031490j,
        -4.28085548650291806830723860758 + 0.96683279082320198340559106474j
    ]

    def assert_expected_values(self, actual):
        self.assert_approx_equal(self.expected_area, actual['area'])
        self.assert_approx_equal(self.expected_imt[0], actual['perimeter'])
        for s in range(2, 7):
            self.assert_approx_equal(self.expected_imt[s], actual['psi%i' % s])
            expected_qs = abs(self.expected_imt[s]) / self.expected_imt[0]
            self.assert_approx_equal(expected_qs, actual['q%i' % s])

    def test_nonconvex_polygon(self):
        minkval = pypaya2.imt_for_polygon(self.vertices)
        self.assert_expected_values(minkval)
        expected_metrics = \
            ['psi%i' % i for i in range(2, 13)] + \
            ['q%i' % i for i in range(2, 13)] + \
            ['area', 'perimeter']
        self.assert_sets_equal(expected_metrics, minkval.keys())

    def test_nonconvex_polygon__from_numpy_array(self):
        minkval = pypaya2.imt_for_polygon(np.array(self.vertices))
        self.assert_expected_values(minkval)
        expected_metrics = \
            ['psi%i' % i for i in range(2, 13)] + \
            ['q%i' % i for i in range(2, 13)] + \
            ['area', 'perimeter']
        self.assert_sets_equal(expected_metrics, minkval.keys())

    def test_missing_argument(self):
        with self.assertRaises(TypeError):
            pypaya2.imt_for_polygon()

    def test_kwargs(self):
        with self.assertRaises(TypeError):
            pypaya2.imt_for_polygon(not_really = 'acceptable')

    def assert_raises_both_ways(self, klass, message, value):
        with self.assertRaises(klass) as cm:
            pypaya2.imt_for_polygon(value)
        self.assertEqual(message, str(cm.exception))
        with self.assertRaises(klass) as cm:
            pypaya2.imt_for_polygon(np.array(value))
        self.assertEqual(message, str(cm.exception))

    def test_bad_argument__empty(self):
        self.assert_raises_both_ways(ValueError, 'data must be 2D, is 0', [[]])

    def test_bad_argument__too_short(self):
        self.assert_raises_both_ways(ValueError, 'polygon must contain at least two vertices', [[0, 0]])

    def test_bad_argument__zero_dimensional(self):
        self.assert_raises_both_ways(ValueError, 'data must be 2D, is 1', [[1.]])

    def test_bad_argument__one_dimensional(self):
        self.assert_raises_both_ways(ValueError, 'data must be 2D, is 1', [[0], [0]])

    def test_bad_argument__three_dimensional(self):
        self.assert_raises_both_ways(ValueError, 'data must be 2D, is 3', [[0, 1, 2], [0, 1, 2]])

if __name__ == "__main__":
    unittest.main()
