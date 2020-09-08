#include "catch.hpp"
#include <papaya2.hpp>
#include <papaya2/pingfile.hpp>
#include "test_helpers.hpp"

using namespace papaya2;

TEST_CASE("width and height can be extracted")
{
    GIVEN("a photo")
    {
        Photo p;
        WHEN("coordinates are set")
        {
            p.set_coordinates(0, 0, 1, 3, 20, 40);
            THEN("size can be extracted")
            {
                CHECK(p.width() == 20);
                CHECK(p.height() == 40);
                vec_t const pix_diag = pixel_diagonal(p);
                CHECK(pix_diag[0] == 1. / 20);
                CHECK(pix_diag[1] == 3. / 40);
            }
        }
    }
}

// scale for TestPhoto
// choose something different from 1 so unit bugs are detectable.
static const double TP_PIXEL_SIDE = 2.;
static const double TP_PIXEL_AREA = fsq(TP_PIXEL_SIDE);

// simple photo struct for use in test cases
struct TestPhoto : BasicPhoto<int>
{
    TestPhoto(const char *data) : TestPhoto(string(data)) {}

    TestPhoto(const string &data)
    {
        // factor image size
        size_t first__ = data.find_first_of('_');
        if (first__ == data.npos)
            throw std::runtime_error("no rowbreak in TestPhoto");
        int numx = first__;
        if (data.size() % (numx + 1u))
            throw std::runtime_error("invalid TestPhoto");
        int numy = data.size() / (numx + 1u);
        if (numx == 0u || numy == 0u)
            throw std::runtime_error("empty TestPhoto");
        set_coordinates(0., 0., TP_PIXEL_SIDE * numx, TP_PIXEL_SIDE * numy,
                        numx, numy);
        // check data
        for (int y = 0; y != numy; ++y) {
            for (int x = 0; x != numx; ++x) {
                int pixel = unchecked_at(data, x, y);
                if (pixel != ' ' /* black */ && pixel != 'x' /* white */)
                    throw std::runtime_error("invalid pixel in TestPhoto");
                else
                    at(x, y) = (pixel == 'x');
            }
            if (unchecked_at(data, numx, y) != '_' /* newline */)
                throw std::runtime_error("missing rowbreak in TestPhoto");
        }
    }

  private:
    int unchecked_at(const string &data, int x, int y)
    {
        return data[(height() - 1 - y) * (width() + 1) + x];
    }
};

TEST_CASE("marching squares algorithm")
{
    // volumes of all the 16 configs at the reference threshold
    // one half
    double const ref_volume[16] = {0.,  .125, .125, .5,  .125, .5,
                                   .75, .875, .125, .75, .5,   .875,
                                   .5,  .875, .875, 1.};
    char const *const config[16] = {"  _  _", "  _x _", "x _  _", "x _x _",
                                    "  _ x_", "  _xx_", "x _ x_", "x _xx_",
                                    " x_  _", " x_x _", "xx_  _", "xx_x _",
                                    " x_ x_", " x_xx_", "xx_ x_", "xx_xx_"};

    SECTION("all 2x2 neighborhoods")
    {
        MinkowskiAccumulator w_imt, b_imt;
        for (int i = 0; i != 16; ++i) {
            INFO("marching squares configuration " << i);
            TestPhoto ph(config[i]);

            SECTION("fundamentals")
            {
                CHECK(ph.width() == 2);
                CHECK(ph.height() == 2);
                CHECK_NOTHROW(check_msq_positively_oriented(ph, CONNECT_WHITE));
                CHECK_NOTHROW(check_msq_positively_oriented(ph, CONNECT_BLACK));
            }
            SECTION("regular marching squares at theshold .75")
            {
                w_imt = imt_regular_marching_squares(ph, .75);
                CHECK(w_imt.area() == approx(ref_volume[i] * TP_PIXEL_AREA));
            }
            SECTION("regular marching squares at threshold 1")
            {
                // connect white case
                b_imt = imt_regular_marching_squares(
                    ph, 1., ANALYZE_BLACK | CONNECT_WHITE);
                w_imt = imt_regular_marching_squares(
                    ph, 1., ANALYZE_WHITE | CONNECT_WHITE);
                CHECK(b_imt.area() + w_imt.area() == approx(TP_PIXEL_AREA));
                CHECK(b_imt.imt(0) == approx(+w_imt.imt(0)));
                CHECK(b_imt.imt(2) == approx(+w_imt.imt(2)));
                CHECK(b_imt.imt(3) == approx(-w_imt.imt(3)));
                // the connect black case
                b_imt = imt_regular_marching_squares(
                    ph, 1., ANALYZE_BLACK | CONNECT_BLACK);
                w_imt = imt_regular_marching_squares(
                    ph, 1., ANALYZE_WHITE | CONNECT_BLACK);
                CHECK(b_imt.area() + w_imt.area() == approx(TP_PIXEL_AREA));
                CHECK(b_imt.imt(0) == approx(+w_imt.imt(0)));
                CHECK(b_imt.imt(2) == approx(+w_imt.imt(2)));
                CHECK(b_imt.imt(3) == approx(-w_imt.imt(3)));
                b_imt = imt_regular_marching_squares(ph, .75);
                CHECK(b_imt.area() == approx(ref_volume[i] * TP_PIXEL_AREA));
            }
            SECTION("interpolated marching squares at threshold .5")
            {
                w_imt = imt_interpolated_marching_squares(ph, .5);
                CHECK(w_imt.area() == approx(ref_volume[i] * TP_PIXEL_AREA));
            }
            SECTION("interpolated marching squares at threshold 1")
            {
                // connect white case
                b_imt = imt_interpolated_marching_squares(
                    ph, 1., ANALYZE_BLACK | CONNECT_WHITE);
                w_imt = imt_interpolated_marching_squares(
                    ph, 1., ANALYZE_WHITE | CONNECT_WHITE);
                CHECK(b_imt.area() + w_imt.area() == approx(TP_PIXEL_AREA));
                CHECK(b_imt.imt(0) == approx(+w_imt.imt(0)));
                CHECK(b_imt.imt(2) == approx(+w_imt.imt(2)));
                CHECK(b_imt.imt(3) == approx(-w_imt.imt(3)));
                // the connect black case
                b_imt = imt_interpolated_marching_squares(
                    ph, 1., ANALYZE_BLACK | CONNECT_BLACK);
                w_imt = imt_interpolated_marching_squares(
                    ph, 1., ANALYZE_WHITE | CONNECT_BLACK);
                CHECK(b_imt.area() + w_imt.area() == approx(TP_PIXEL_AREA));
                CHECK(b_imt.imt(0) == approx(+w_imt.imt(0)));
                CHECK(b_imt.imt(2) == approx(+w_imt.imt(2)));
                CHECK(b_imt.imt(3) == approx(-w_imt.imt(3)));
            }
        }
    }

    SECTION("negating should flip the normals")
    {
        Photo ph;
        ph.set_coordinates(0, 0, 2, 2, 2, 2);
        ph(0, 1) = -1.3;
        ph(1, 1) = 0.7;
        ph(0, 0) = -2.1;
        ph(1, 0) = 1.7;
        SECTION("config #12")
        {
            auto i_imt = imt_interpolated_marching_squares(ph, 0.);
            auto b_imt =
                imt_interpolated_marching_squares(ph, 0., ANALYZE_BLACK);
            CHECK(i_imt.imt(3) == approx(-b_imt.imt(3)));
            CHECK(i_imt.imt(0) == approx(b_imt.imt(0)));
        }
        SECTION("config #12, connect black")
        {
            auto i_imt =
                imt_interpolated_marching_squares(ph, 0., CONNECT_BLACK);
            auto b_imt = imt_interpolated_marching_squares(
                ph, 0., ANALYZE_BLACK | CONNECT_BLACK);
            CHECK(i_imt.imt(3) == approx(-b_imt.imt(3)));
            CHECK(i_imt.imt(0) == approx(b_imt.imt(0)));
        }
        ph(0, 1) = 0.7;
        ph(1, 1) = -1.3;
        ph(0, 0) = -2.1;
        ph(1, 0) = 1.7;
        SECTION("critical config #6")
        {
            auto i_imt = imt_interpolated_marching_squares(ph, 0.);
            auto b_imt =
                imt_interpolated_marching_squares(ph, 0., ANALYZE_BLACK);
            CHECK(i_imt.imt(3) == approx(-b_imt.imt(3)));
            CHECK(i_imt.imt(0) == approx(b_imt.imt(0)));
        }
        SECTION("critical config #6, connect white")
        {
            auto i_imt =
                imt_interpolated_marching_squares(ph, 0., CONNECT_BLACK);
            auto b_imt = imt_interpolated_marching_squares(
                ph, 0., ANALYZE_BLACK | CONNECT_BLACK);
            CHECK(i_imt.imt(3) == approx(-b_imt.imt(3)));
            CHECK(i_imt.imt(0) == approx(b_imt.imt(0)));
        }
        ph(0, 1) = -0.7;
        ph(1, 1) = 1.3;
        ph(0, 0) = 2.1;
        ph(1, 0) = -1.7;
        SECTION("critical config #9")
        {
            auto i_imt = imt_interpolated_marching_squares(ph, 0.);
            auto b_imt =
                imt_interpolated_marching_squares(ph, 0., ANALYZE_BLACK);
            CHECK(i_imt.imt(3) == approx(-b_imt.imt(3)));
            CHECK(i_imt.imt(0) == approx(b_imt.imt(0)));
        }
        SECTION("critical config #9, connect white")
        {
            auto i_imt =
                imt_interpolated_marching_squares(ph, 0., CONNECT_BLACK);
            auto b_imt = imt_interpolated_marching_squares(
                ph, 0., ANALYZE_BLACK | CONNECT_BLACK);
            CHECK(i_imt.imt(3) == approx(-b_imt.imt(3)));
            CHECK(i_imt.imt(0) == approx(b_imt.imt(0)));
        }
    }
    SECTION("a single pixel")
    {
        TestPhoto one_pixel = "   _"
                              " x _"
                              "   _";
        double const ref_area = (.5) * TP_PIXEL_AREA;
        double const ref_peri = (2 * SQRT2) * TP_PIXEL_SIDE;
        double const ref_psi2 = (0.) * TP_PIXEL_SIDE;
        double const ref_psi3 = (0.) * TP_PIXEL_SIDE;
        double const ref_psi4 = (-2 * SQRT2) * TP_PIXEL_SIDE;
        auto i_imt = imt_interpolated_marching_squares(one_pixel, .5);
        auto r_imt = imt_regular_marching_squares(one_pixel, .75);
        CHECK(i_imt.area() == approx(ref_area));
        CHECK(i_imt.msm(2) == approx(0.));
        CHECK(i_imt.msm(3) == approx(0.));
        CHECK(i_imt.msm(4) == approx(1.));
        CHECK(i_imt.imt(0) == approx(ref_peri));
        CHECK(i_imt.imt(2) == approx(ref_psi2));
        CHECK(i_imt.imt(3) == approx(ref_psi3));
        CHECK(i_imt.imt(4) == approx(ref_psi4));
        CHECK(r_imt.area() == approx(ref_area));
        CHECK(r_imt.msm(2) == approx(0.));
        CHECK(r_imt.msm(3) == approx(0.));
        CHECK(r_imt.msm(4) == approx(1.));
        CHECK(r_imt.imt(0) == approx(ref_peri));
        CHECK(r_imt.imt(2) == approx(ref_psi2));
        CHECK(r_imt.imt(3) == approx(ref_psi3));
        CHECK(r_imt.imt(4) == approx(ref_psi4));
        r_imt = imt_regular_marching_squares(one_pixel, 1);
        CHECK(r_imt.area() == approx(ref_area));
        r_imt = imt_regular_marching_squares(one_pixel, 1, ANALYZE_BLACK);
        CHECK(r_imt.area() == approx(4 * (1. - .125) * TP_PIXEL_AREA));
    }
    SECTION("tetromino1")
    {
        // reference values in pixels
        auto const ref_area = 3.5 * TP_PIXEL_AREA;
        auto const ref_peri = (4 * SQRT2 + 2) * TP_PIXEL_SIDE;

        SECTION("orientation #1")
        {
            auto const ref_psi2 = complex_t(2., 0.) * TP_PIXEL_SIDE;
            auto const ref_psi3 = complex_t(-4., 0.) * TP_PIXEL_SIDE;
            auto const ref_psi4 = complex_t(2 - 4 * SQRT2, 0.) * TP_PIXEL_SIDE;
            TestPhoto ph = "    _"
                           " x  _"
                           " xx _"
                           " x  _"
                           "    _";
            auto i_imt = imt_interpolated_marching_squares(ph, .5);
            auto r_imt = imt_regular_marching_squares(ph, .75);
            CHECK(i_imt.area() == approx(ref_area));
            CHECK(i_imt.perimeter() == approx(ref_peri));
            CHECK(i_imt.imt(0) == approx(ref_peri));
            CHECK(i_imt.imt(2) == approx(ref_psi2));
            CHECK(i_imt.imt(3) == approx(ref_psi3));
            CHECK(i_imt.imt(4) == approx(ref_psi4));
            CHECK(r_imt.area() == approx(ref_area));
            CHECK(r_imt.perimeter() == approx(ref_peri));
            CHECK(r_imt.imt(0) == approx(ref_peri));
            CHECK(r_imt.imt(2) == approx(ref_psi2));
            CHECK(r_imt.imt(3) == approx(ref_psi3));
            CHECK(r_imt.imt(4) == approx(ref_psi4));
            i_imt = imt_interpolated_marching_squares(ph, .5, ANALYZE_BLACK);
            r_imt = imt_regular_marching_squares(ph, .75, ANALYZE_BLACK);
            CHECK(i_imt.area() == approx(12 * TP_PIXEL_AREA - ref_area));
            CHECK(i_imt.perimeter() == approx(ref_peri));
            CHECK(i_imt.imt(0) == approx(ref_peri));
            CHECK(i_imt.imt(2) == approx(ref_psi2));
            CHECK(i_imt.imt(3) == approx(-ref_psi3));
            CHECK(i_imt.imt(4) == approx(ref_psi4));
            CHECK(r_imt.imt(2) == approx(ref_psi2));
            CHECK(r_imt.imt(3) == approx(-ref_psi3));
            CHECK(r_imt.imt(4) == approx(ref_psi4));
        }

        SECTION("orientation #2")
        {
            auto const ref_psi2 = complex_t(-2., 0.) * TP_PIXEL_SIDE;
            auto const ref_psi3 = complex_t(0., 4.) * TP_PIXEL_SIDE;
            auto const ref_psi4 = complex_t(2 - 4 * SQRT2, 0.) * TP_PIXEL_SIDE;
            SECTION("test w/o padding")
            {

                TestPhoto ph = "     _"
                               "  x  _"
                               " xxx _"
                               "     _";
                auto i_imt = imt_interpolated_marching_squares(ph, .5);
                auto r_imt = imt_regular_marching_squares(ph, .75);
                CHECK(i_imt.area() == approx(ref_area));
                CHECK(i_imt.perimeter() == approx(ref_peri));
                CHECK(i_imt.imt(2) == approx(ref_psi2));
                CHECK(i_imt.imt(3) == approx(ref_psi3));
                CHECK(i_imt.imt(4) == approx(ref_psi4));
                CHECK(r_imt.area() == approx(ref_area));
                CHECK(r_imt.perimeter() == approx(ref_peri));
                CHECK(r_imt.imt(2) == approx(ref_psi2));
                CHECK(r_imt.imt(3) == approx(ref_psi3));
                CHECK(r_imt.imt(4) == approx(ref_psi4));
                i_imt =
                    imt_interpolated_marching_squares(ph, .5, ANALYZE_BLACK);
                r_imt = imt_regular_marching_squares(ph, .75, ANALYZE_BLACK);
                CHECK(i_imt.area() == approx(12 * TP_PIXEL_AREA - ref_area));
                CHECK(i_imt.perimeter() == approx(ref_peri));
                CHECK(i_imt.imt(0) == approx(ref_peri));
                CHECK(i_imt.imt(2) == approx(ref_psi2));
                CHECK(i_imt.imt(3) == approx(-ref_psi3));
                CHECK(i_imt.imt(4) == approx(ref_psi4));
                CHECK(r_imt.imt(2) == approx(ref_psi2));
                CHECK(r_imt.imt(3) == approx(-ref_psi3));
                CHECK(r_imt.imt(4) == approx(ref_psi4));
            }
            SECTION("orientation #2, test padding")
            {
                TestPhoto te = "    _"
                               " x  _"
                               "xxx _";
                auto const ph = make_padded_view(te, 0.);
                auto const i_imt = imt_interpolated_marching_squares(ph, .5);
                auto const r_imt = imt_regular_marching_squares(ph, .75);
                CHECK(i_imt.area() == approx(ref_area));
                CHECK(i_imt.perimeter() == approx(ref_peri));
                CHECK(i_imt.imt(2) == approx(ref_psi2));
                CHECK(i_imt.imt(3) == approx(ref_psi3));
                CHECK(i_imt.imt(4) == approx(ref_psi4));
                CHECK(r_imt.area() == approx(ref_area));
                CHECK(r_imt.perimeter() == approx(ref_peri));
                CHECK(r_imt.imt(2) == approx(ref_psi2));
                CHECK(r_imt.imt(3) == approx(ref_psi3));
                CHECK(r_imt.imt(4) == approx(ref_psi4));
            }
            SECTION("orientation #2, test other padding")
            {
                TestPhoto te = "    _"
                               "  x _"
                               " xxx_";
                auto const ph = make_padded_view(te, 0.);
                auto const i_imt = imt_interpolated_marching_squares(ph, .5);
                auto const r_imt = imt_regular_marching_squares(ph, .75);
                CHECK(i_imt.area() == approx(ref_area));
                CHECK(i_imt.perimeter() == approx(ref_peri));
                CHECK(i_imt.imt(2) == approx(ref_psi2));
                CHECK(i_imt.imt(3) == approx(ref_psi3));
                CHECK(i_imt.imt(4) == approx(ref_psi4));
                CHECK(r_imt.area() == approx(ref_area));
                CHECK(r_imt.perimeter() == approx(ref_peri));
                CHECK(r_imt.imt(2) == approx(ref_psi2));
                CHECK(r_imt.imt(3) == approx(ref_psi3));
                CHECK(r_imt.imt(4) == approx(ref_psi4));
            }
        }
    }
    SECTION("corner case where the threshold matches a pixel value exactly")
    {
        TestPhoto one_pixel = "   _"
                              " x _"
                              "   _";
        auto i_imt = imt_interpolated_marching_squares(one_pixel, 1.);
        CHECK(i_imt.area() == 0.);
        CHECK(i_imt.imt(0) == 0.);
        CHECK(i_imt.imt(2) == 0.);
        CHECK(i_imt.imt(3) == 0.);
        CHECK(i_imt.imt(4) == 0.);
        i_imt = imt_interpolated_marching_squares(one_pixel, 1., ANALYZE_BLACK);
        CHECK(i_imt.area() == 4. * TP_PIXEL_AREA);
        CHECK(i_imt.imt(0) == 0.);
        CHECK(i_imt.imt(2) == 0.);
        CHECK(i_imt.imt(3) == 0.);
        CHECK(i_imt.imt(4) == 0.);
    }
    SECTION("Minkowski maps")
    {
        TestPhoto one_pixel = "x_";
        auto const padded = make_padded_view(one_pixel, 0.);

        complex_image_t out;
        minkowski_map_interpolated_marching_squares(&out, padded, .5, 2);
        CHECK(out.width() == 2);
        CHECK(out.height() == 2);
        CHECK(out(0, 0) != 0.);
        CHECK(out(0, 1) != 0.);
        CHECK(out(1, 1) != 0.);
        CHECK(out(1, 0) != 0.);
        minkowski_map_interpolated_marching_squares(&out, one_pixel, .5, 2);
        CHECK(out.width() == 0);
        CHECK(out.height() == 0);
    }
    SECTION("Grayscale image")
    {
        PingFile potato("validation_data/potato.png");
        auto const imt =
            imt_interpolated_marching_squares(potato, 1.469734492275599e+02);
        using std::arg;
        CHECK(imt.area() == approx(7.760018008530601e+04));
        CHECK(imt.perimeter() == approx(1.009820521813200e+03));
        CHECK(imt.msm(2) == approx(1.563708314579508e-01));
        CHECK(arg(imt.imt(2)) == approx(-5.574818755549115e-1));
        CHECK(imt.msm(3) == approx(3.980635759575531e-02));
        CHECK(arg(imt.imt(3)) == approx(4.577327494757427e-01));
        CHECK(imt.msm(4) == approx(1.481665040254243e-01));
        CHECK(arg(imt.imt(4)) == approx(7.141179567742937e-01));
        CHECK(imt.msm(5) == approx(1.922804813953316e-01));
        CHECK(arg(imt.imt(5)) == approx(2.070654574313143e+00));
        CHECK(imt.msm(6) == approx(1.454253098691293e-01));
        CHECK(arg(imt.imt(6)) == approx(-2.519436436486195e-1));
        CHECK(imt.msm(7) == approx(7.671564980492507e-03));
        CHECK(arg(imt.imt(7)) == approx(1.404644402590062e+00));
        CHECK(imt.msm(8) == approx(1.860509319491248e-01));
        CHECK(arg(imt.imt(8)) == approx(7.356981574804343e-02));
    }
}
