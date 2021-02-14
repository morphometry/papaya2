#!/usr/bin/env python3
from test_helper import *
import numpy as np
from PIL import Image

class ImageTest(PypayaTestCase):
    def test_imt_for_image__computes_accurately(self):
        image = self.image_fixture('potato.png')
        self.assert_equal((302, 361), image.shape)

        minkval = pypaya2.imt_for_image(image, threshold = 1.469734492275599e+02)

        self.expected_area = 7.760018008530601e+04
        def imt_from_polar(q_s, arg):
            return 1.009820521813200e+03 * q_s * (np.cos(arg) + 1j * np.sin(arg))
        self.expected_imt = [
            imt_from_polar(1, 0),
            imt_from_polar(0, 0),
            imt_from_polar(1.563708314579508e-01, -5.574818755549115e-1),
            imt_from_polar(3.980635759575531e-02, 4.577327494757427e-01),
            imt_from_polar(1.481665040254243e-01, 7.141179567742937e-01),
            imt_from_polar(1.922804813953316e-01, 2.070654574313143e+00),
            imt_from_polar(1.454253098691293e-01, -2.519436436486195e-1),
            imt_from_polar(7.671564980492507e-03, 1.404644402590062e+00),
            imt_from_polar(1.860509319491248e-01, 7.356981574804343e-02)
        ]

        self.assert_expected_values(minkval, accuracy = 1e-12)

    def test_imt_for_image__without_threshold__defaults_to_one_half(self):
        image = self.image_fixture('potato.png')
        self.assert_equal(0., image.min())
        self.assert_equal(255., image.max())
        image /= 255.

        minkval = pypaya2.imt_for_image(image)

        self.expected_area = 7.76960185776013240684e+04
        self.expected_imt = np.array([
            1.00920193792553550338e+03,
            0.,
            1.33895776388211203312e+02 - 8.41864463565138834156e+01j,
            3.67446495067691216718e+01 + 1.94915084343996412031e+01j,
            1.10618461907130892996e+02 + 1.00873659621835642497e+02j,
            -9.5708454838491164196e+01 + 1.78202907585352107844e+02j,
            1.49139226616021858263e+02 - 4.23062941599343318444e+01j,
            1.11609328312643918402e-01 + 3.12616420762001157030e+00j,
            1.58742464698595824757e+02 + 1.43409448675481456803e+01j
        ])

        self.assert_expected_values(minkval, accuracy = 1e-12)

    def test_imt_for_image__missing_argument(self):
        with self.assertRaises(TypeError):
            pypaya2.imt_for_image()

    def test_kwargs(self):
        with self.assertRaises(ValueError):
            pypaya2.imt_for_image(self.mock_image(), not_really = 'acceptable')

    def test_imt_for_image__bad_argument__vector(self):
        self.assert_raises_both_ways(ValueError, '# dimensions of image must be 2, is 1', [0, 1])

    def test_imt_for_image__bad_argument__tensor(self):
        self.assert_raises_both_ways(ValueError, '# dimensions of image must be 2, is 3', [[[0, 1], [0, 1]], [[0, 1], [0, 1]]])

    def test_pil__image_orientation(self):
        coords_test = self.image_fixture('coordinates.png')
        self.assert_equal((2, 3), coords_test.shape)
        self.assert_equal(255., coords_test[0,0])
        self.assert_equal(255., coords_test[1,0])
        self.assert_equal(120., coords_test[0,1])
        self.assert_equal(240., coords_test[1,1])
        self.assert_equal(0.,   coords_test[0,2])
        self.assert_equal(60.,  coords_test[1,2])

    def test_pil__image_orientation_and_channels(self):
        coords_test = self.image_fixture('channels.png')
        self.assert_equal((2, 3), coords_test.shape)
        self.assert_equal(255., coords_test[0,0])
        self.assert_equal(255., coords_test[1,0])
        self.assert_equal(120., coords_test[0,1])
        self.assert_equal(240., coords_test[1,1])
        self.assert_equal(0.,   coords_test[0,2])
        self.assert_equal(60.,  coords_test[1,2])

    def mock_image(self):
        return np.eye(2)

    def image_fixture(self, name):
        pil_image = Image.open('../../test/validation_data/' + name)
        image = np.asarray(pil_image)
        if len(image.shape) == 3:
            # take red channel
            image = image[:,:,0]
        # swap coordinates s.t. image(x,y)
        image = np.transpose(image)
        # flip y axis
        image = image[:,-1::-1]
        return image.astype(np.float64)

    def assert_expected_values(self, actual, accuracy = 1e-14):
        self.assert_approx_equal(self.expected_area, actual['area'], accuracy = accuracy)
        self.assert_approx_equal(self.expected_imt[0], actual['perimeter'], accuracy = accuracy)
        for s in range(2, len(self.expected_imt)):
            self.assert_approx_equal(self.expected_imt[s], actual['psi%i' % s], accuracy = accuracy)
            expected_qs = abs(self.expected_imt[s]) / self.expected_imt[0]
            self.assert_approx_equal(expected_qs, actual['q%i' % s], accuracy = accuracy)

    def assert_raises_both_ways(self, klass, message, value):
        with self.assertRaises(klass) as cm:
            pypaya2.imt_for_image(value)
        self.assertEqual(message, str(cm.exception))
        with self.assertRaises(klass) as cm:
            pypaya2.imt_for_image(np.asarray(value))
        self.assertEqual(message, str(cm.exception))

if __name__ == "__main__":
    unittest.main()
