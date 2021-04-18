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

    def test_imt_for_image__multiple_thresholds(self):
        image = self.image_fixture('potato.png')
        minkval = pypaya2.imt_for_image(image, threshold = [1.469734492275599e+02, 127.5])

        accuracy = 1e-12
        self.assert_approx_equal(7.760018008530601e+04, minkval['area'][0], accuracy = accuracy)
        self.assert_approx_equal(7.769601857760132e+04, minkval['area'][1], accuracy = accuracy)

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

    def test_imt_for_image__bad_argument__tensor_threshold(self):
        with self.assertRaises(ValueError) as cm:
            pypaya2.imt_for_image(self.mock_image(), threshold = np.eye(3))
        self.assert_equal('cannot cast higher-dimensional array to vector for threshold argument', str(cm.exception))

    def test_imt_for_image__bad_argument__string_threshold(self):
        with self.assertRaises(ValueError) as cm:
            pypaya2.imt_for_image(self.mock_image(), threshold = 'xyz')
        self.assert_equal('cannot cast data for threshold argument', str(cm.exception))

    def test_minkowski_map_for_image(self):
        image = np.zeros(shape=(2, 3))
        image[:, 1] = 1.
        minkmap = pypaya2.minkowski_map_for_image(image, boundary='periodic')
        self.assert_equal((1, 2, 3), minkmap.shape)
        expected = [[-1, -1, 0]] * 2
        self.assert_approx_equal(expected, minkmap[0], accuracy=1e-6)

    def test_minkowski_map_for_image_different_s(self):
        image = np.zeros(shape=(2, 3))
        image[:, 1] = 1.
        with self.assertRaises(ValueError) as cm:
          minkmap = pypaya2.minkowski_map_for_image(image, boundary='periodic', s = [0,2])
        self.assert_equal('illegal value for s argument, must be a single int value', str(cm.exception))

        with self.assertRaises(ValueError) as cm:
          minkmap = pypaya2.minkowski_map_for_image(image, boundary='periodic', s = 0.67)
        self.assert_equal('illegal value for s argument, s has to be an int', str(cm.exception))

        with self.assertRaises(ValueError) as cm:
          minkmap = pypaya2.minkowski_map_for_image(image, boundary='periodic', s = -1)
        self.assert_equal('illegal value for s argument, s has to be a non-negative int', str(cm.exception))

        minkmap = pypaya2.minkowski_map_for_image(image, boundary='periodic', s = 0)
        self.assert_equal((1, 2, 3), minkmap.shape)
        expected = [[1, 1, 0]] * 2
        self.assert_approx_equal(expected, minkmap[0], accuracy=1e-6)

        with self.assertRaises(ValueError) as cm:
          minkmap = pypaya2.minkowski_map_for_image(image, boundary='periodic', s = 1)
        self.assert_equal('illegal value for s argument, s cannot be 1', str(cm.exception))

        minkmap = pypaya2.minkowski_map_for_image(image, boundary='periodic', s = 2)
        self.assert_equal((1, 2, 3), minkmap.shape)
        expected = [[-1, -1, 0]] * 2
        self.assert_approx_equal(expected, minkmap[0], accuracy=1e-6)

        minkmap = pypaya2.minkowski_map_for_image(image, boundary='periodic', s = 3)
        self.assert_equal((1, 2, 3), minkmap.shape)
        expected = [[0+1j, 0-1j, 0]] * 2
        self.assert_approx_equal(expected, minkmap[0], accuracy=1e-6)

        with self.assertRaises(ValueError) as cm:
          minkmap = pypaya2.minkowski_map_for_image(image, boundary='periodic', s = 67)
        self.assert_equal('illegal value for s argument, s needs to be smaller than MAX_S', str(cm.exception))

    def test_minkowski_map_for_image__padded__zeros(self):
        image = self.mock_image()
        self.assert_equal((3, 2), image.shape)
        minkmap = pypaya2.minkowski_map_for_image(image, threshold=0.06)
        self.assert_equal((1, 4, 3), minkmap.shape)
        expected = [
            [+0.00000000e+00+0.98994949j, +9.32903324e-01+0.41905818j, -2.87102621e-16-1.29299526j],
            [-9.66685285e-01+0.29668091j, +0.00000000e+00+0.j,         -9.80555919e-01+0.22709318j],
            [-7.27257400e-01-0.82072935j, +0.00000000e+00+0.j,         -9.85086818e-01-0.19900744j],
            [-3.45420341e-16-0.56568542j, +6.70820393e-01-0.89442719j, +0.00000000e+00+1.27279221j]
        ]
        self.assert_approx_equal(expected, minkmap[0], accuracy=1e-6)

    def test_minkowski_map_for_image__padded__padding_value(self):
        image = self.mock_image()
        self.assert_equal((3, 2), image.shape)
        minkmap = pypaya2.minkowski_map_for_image(image, threshold=0.06, boundary=0.03)
        self.assert_equal((1, 4, 3), minkmap.shape)
        expected = [
            [+0.00000000e+00+1.16464646j, 9.74244511e-01+0.26113419j, +2.99957962e-16-1.35089057j],
            [-9.86423194e-01+0.18991693j, 0.00000000e+00+0.j,         -9.93416388e-01+0.13237905j],
            [-8.30535726e-01-0.6564754j,  0.00000000e+00+0.j,         -9.94880423e-01-0.11675958j],
            [-1.79439138e-16-0.80812204j, 8.03748429e-01-0.7037892j,  +0.00000000e+00+1.33978127j]
        ]
        self.assert_approx_equal(expected, minkmap[0], accuracy=1e-6)

    def test_minkowski_map_for_image__periodic(self):
        image = self.mock_image()
        self.assert_equal((3, 2), image.shape)
        minkmap = pypaya2.minkowski_map_for_image(image, threshold=0.6, boundary='periodic')
        self.assert_equal((1, 3, 2), minkmap.shape)
        expected = [
            [-0.07027819-0.31234752j, -0.07027819+0.31234752j],
            [+0.00000000+0.j,         +0.00000000+0.j],
            [-0.94135745+0.39223227j, -0.94135745-0.39223227j]
        ]
        self.assert_approx_equal(expected, minkmap[0], accuracy=1e-6)

    def test_minkowski_map_for_image__multi_threshold(self):
        image = self.mock_image()
        self.assert_equal((3, 2), image.shape)
        minkmap = pypaya2.minkowski_map_for_image(image, threshold=[0.6, 1.1], boundary='periodic')
        self.assert_equal((2, 3, 2), minkmap.shape)
        expected = [
            [-0.07027819-0.31234752j, -0.07027819+0.31234752j],
            [+0.00000000+0.j,         +0.00000000+0.j],
            [-0.94135745+0.39223227j, -0.94135745-0.39223227j]
        ]
        self.assert_approx_equal(expected, minkmap[0], accuracy=1e-6)
        expected = np.zeros(shape=(3, 2))
        self.assert_approx_equal(expected, minkmap[1], accuracy=1e-6)

    def test_minkowski_map_for_image__bad_threshold_argument(self):
        with self.assertRaises(ValueError) as cm:
            pypaya2.minkowski_map_for_image(self.mock_image(), threshold='nonsense')
        self.assert_equal('cannot cast data for threshold argument', str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            pypaya2.minkowski_map_for_image(self.mock_image(), threshold=ValueError())
        self.assert_equal('cannot cast data for threshold argument', str(cm.exception))

    def test_minkowski_map_for_image__bad_boundary_argument(self):
        with self.assertRaises(ValueError) as cm:
            pypaya2.minkowski_map_for_image(self.mock_image(), boundary='nonsense')
        self.assert_equal('illegal value for boundary keyword argument: nonsense', str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            pypaya2.minkowski_map_for_image(self.mock_image(), boundary=ValueError())
        self.assert_equal('cannot cast data for boundary argument', str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            pypaya2.minkowski_map_for_image(self.mock_image(), boundary=[1., 2.])
        self.assert_equal("illegal value for boundary keyword argument, must be a Float or the string 'periodic'", str(cm.exception))

    def test_minkowski_map_for_image__bad_keyword_argument(self):
        with self.assertRaises(ValueError) as cm:
            pypaya2.minkowski_map_for_image(self.mock_image(), no_such='argument')
        self.assert_equal('illegal keyword argument: no_such', str(cm.exception))

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
        return np.transpose([
            [0.2, 0.4, 0.1],
            [0.7, 0.3, 0.6],
        ])

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
