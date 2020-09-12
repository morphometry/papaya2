#!/usr/bin/env python3
from test_helper import *
import numpy as np

class PointPatternTest(PypayaTestCase):
    def test_basic_point_pattern_analysis(self):
        x = np.random.uniform(size = [100, 2])
        minkval = pypaya2.imt_for_pointpattern(x)
        expected = ['area', 'perimeter', 'q2', 'psi2', 'q3', 'psi3', 'q4', 'psi4', 'q5', 'psi5', 'q6', 'psi6', 'q7', 'psi7', 'q8', 'psi8', 'q9', 'psi9', 'q10', 'psi10', 'q11', 'psi11', 'q12', 'psi12']
        self.assert_sets_equal(expected, (minkval.keys()))

if __name__ == "__main__":
    unittest.main()
