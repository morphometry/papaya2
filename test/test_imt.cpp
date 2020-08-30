#include "catch.hpp"
#include <papaya2.hpp>
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

    SECTION("imt_polyon with concave polygon")
    {
        auto const vertices = std::vector<vec_t>{
            {2., 3.}, {.25, .5}, {-3., 4.}, {-2., 0.}, {2., 0.}};
        auto const expected_area = 67. / 8;
        auto const expected_imt = std::vector<complex_t>{
            complex_t(18.9509878231867863325879476083, 0.),
            complex_t(0., 0.),
            complex_t(4.03586622195516544754949115540, 3.83613097404496685735871370308),
            complex_t(-2.25242470620683687651763840592, 7.91410261263648605530287871331),
            complex_t(2.23638454508277526919840618339, 2.16586100869600377036856029938),
            complex_t(0.70071569649771790984433739720, -10.30885074737037023789259031490),
            complex_t(-4.28085548650291806830723860758, 0.96683279082320198340559106474)
        };

        SECTION("concave polygon, relabel vertices")
        {
            const size_t N = vertices.size();
            std::vector<vec_t> relabelled_vertices(N);
            for (size_t r = 0; r != N; ++r)
            {
                for (size_t i = 0; i != N; ++i)
                    relabelled_vertices.at((i + r) % N) = vertices.at(i);

                auto imt = imt_polygon(relabelled_vertices);
                CHECK(imt.area() == approx(expected_area));
                for (int s : {0, 2, 3, 4, 5, 6})
                    CHECK(imt.imt(s) == approx(expected_imt[s]));
            }
        }

        SECTION("concave polygon, zero-length edges")
        {
            const size_t N = vertices.size();
            std::vector<vec_t> relabelled_vertices(N + 1);
            for (size_t r = 0; r != N; ++r)
            {
                size_t i = 0;
                for (; i != r; ++i)
                    relabelled_vertices.at(i) = vertices.at(i);
                relabelled_vertices.at(i) = vertices.at(r);
                for (; i != N; ++i)
                    relabelled_vertices.at(i + 1) = vertices.at(i);

                auto imt = imt_polygon(relabelled_vertices);
                CHECK(imt.area() == approx(expected_area));
                for (int s : {0, 2, 3, 4, 5, 6})
                    CHECK(imt.imt(s) == approx(expected_imt[s]));
            }
        }
    }

    SECTION("imt_polygon with equilateral triangle")
    {
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
