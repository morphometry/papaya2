#include "catch.hpp"
#include "papaya2.hpp"
#include "test_helpers.hpp"
#include <sstream>

using namespace papaya2;

TEST_CASE("GnuplottableContour struct")
{
    std::ostringstream os;
    GnuplottableContour cont(os);
    add_contour_segment(&cont, {0.25, 0.75}, {0.25, 0.5}, {1., 1.});
    CHECK(os.str() == "5.000000000000000e-01 1.250000000000000e+00 "
                      "7.500000000000000e-01 5.000000000000000e-01 \n");
}

static vec_t ellipse_polar_form(double e, double phi)
{
    const double r = std::pow(1 - fsq(e * std::cos(phi)), -.5);
    return r * cos_sin(phi);
}

static complex_t expi(double phi) { return std::polar(1., phi); }

static void trace_ellipse(MinkowskiAccumulator *acc, double e,
                          double num_segments)
{
    vec_t last, now = ellipse_polar_form(e, 0);
    const int n = num_segments;
    const vec_t off = {0., 0.};
    for (int i = 1; i != n + 1; ++i) {
        last = now;
        now = ellipse_polar_form(e, i * TWO_PI / n);
        add_contour_segment(acc, off, last, now);
    }
}

TEST_CASE("MinkowskiAccumulator struct")
{
    const int num_segments = 30000;
    SECTION("shape is a unit circle")
    {
        MinkowskiAccumulator acc;
        trace_ellipse(&acc, 0., num_segments);
        CHECK(acc.imt(2) == approx(0.));
        CHECK(acc.imt(3) == approx(0.));
        CHECK(acc.imt(4) == approx(0.));
    }
    SECTION("shape is a 3:8 ellipse")
    {
        MinkowskiAccumulator acc;
        trace_ellipse(&acc, std::sqrt(1 - 9. / 64), num_segments);
        CHECK(acc.imt(2) / acc.perimeter() ==
              approx(-0.63070026151848068207).tolerance(1e-6));
        CHECK(acc.imt(3) / acc.perimeter() == approx(0.).tolerance(1e-6));
        CHECK(acc.imt(4) / acc.perimeter() ==
              approx(0.34844502478902471207).tolerance(1e-6));
        CHECK(acc.imt(6) / acc.perimeter() ==
              approx(-0.18211710641944390473).tolerance(1e-6));
        CHECK(acc.imt(8) / acc.perimeter() ==
              approx(0.092302729744212168710).tolerance(1e-6));
    }

    SECTION("large square should approach q4=1")
    {
        const int resolution = 200;
        const double ref_area_pix =
            fsq(resolution) - 4 * (resolution - 1.) - .5;
        const double ref_peri_pix = 4 * (resolution - 3.) + 4 / SQRT2;
        const double ref_area = ref_area_pix * fsq(2. / resolution);
        const double ref_peri = ref_peri_pix * (2. / resolution);
        struct PixelizedSquare : Photo
        {
            PixelizedSquare()
            {
                set_coordinates(-1, -1, 1, 1, resolution, resolution);
                for (int i = 1; i < resolution - 1; ++i) {
                    for (int j = 1; j < resolution - 1; ++j) {
                        at(i, j) = 1;
                    }
                }
            }
        } ph;

        auto imt = imt_interpolated_marching_squares(ph, .5);
        // FIXME accuracy problem?
        CHECK(imt.area() == approx(ref_area).tolerance(1e-6));
        CHECK(imt.perimeter() == approx(ref_peri));
        CHECK(imt.msm(4) > .99);
    }

    SECTION("imt_polygon")
    {
        SECTION("concave polygon")
        {
            auto const vertices = std::vector<vec_t>{
                {2., 3.}, {0., 0.}, {-2., 3.}, {-2., 0.}, {2., 0.}};
            auto imt = imt_polygon(vertices);
            CHECK(imt.area() == approx(6.));
        }

        SECTION("equilateral triangle")
        {
            vec_t const v0 = {1., 0.}, v1 = cos_sin(TWO_PI / 3),
                        v2 = cos_sin(-TWO_PI / 3);
            auto const vertices = std::vector<vec_t>{v0, v1, v2};
            MinkowskiAccumulator imt_tri, imt_poly;
            imt_poly = imt_polygon(vertices);
            add_triangle_area(&imt_tri, ORIGIN, v0, v1, v2);
            // polygon and triangle formulas should be equivalent
            CHECK(imt_poly.area() == approx(imt_tri.area()));
            // check IMT against exact values for equilateral triangle
            double const side = imt_poly.perimeter() / 3;
            for (int s : {2, 3, 4, 5, 6}) {
                complex_t const ref_psi =
                    side * (expi(TWO_PI / 6 * s) + expi(-TWO_PI / 6 * s) +
                            expi(TWO_PI / 2 * s));
                CHECK(imt_poly.imt(s) == approx(ref_psi));
            }
        }
    }
}
