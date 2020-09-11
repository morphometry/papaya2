// sersic: small program which samples an elliptic Sérsic profile and computes
// IMT's for a number of interpolated marching squares thresholds
// using the preliminary Papaya2 library.
// 2018 Sebastian Kapfer <sebastian.kapfer@fau.de>
#include "sersic.hpp"
#include "common.hpp"

using namespace papaya2;

int main(int, const char **argv)
{
    if (!argv[1])
        die("give mode");

    double aspect = 3. / 8;
    void (Photo::*method)(const Sersic &) = &Photo::sample_function;
    Photo ds;
    bool interpolated_marching_squares = false;
    double threshold = 1.2;
    int resolution = 200;

    string mode = argv[1];
    for (argv += 2; *argv; ++argv) {
        if (string(*argv) == "integrate") {
            method = &Photo::integrate_function;
        } else if (string(*argv) == "threshold") {
            threshold = read_arg<double>(argv++);
        } else if (string(*argv) == "aspect") {
            aspect = read_arg<double>(argv++);
        } else if (string(*argv) == "interpolated_marching_squares") {
            interpolated_marching_squares = true;
        } else if (string(*argv) == "resolution") {
            resolution = read_arg<unsigned long>(argv++);
        } else {
            std::cerr << "illegal argument: " << *argv << "\n";
            return 1;
        }
    }

    ds.set_coordinates(-1, -1, 1, 1, resolution, resolution);
    MinkowskiAccumulator imt;

    if (mode == "scan_threshold") {
        (ds.*method)(Sersic(1., aspect, deg2rad(10)));
        Datafile df(std::cout << "# threshold area perimeter q2 q3 q4\n");
        for (auto threshold : logspace(1.1, 2.71828, 30, true)) {
            if (!interpolated_marching_squares)
                imt = imt_regular_marching_squares(ds, threshold);
            else
                imt = imt_interpolated_marching_squares(ds, threshold);
            df << threshold << imt.area() << imt.perimeter() << imt.msm(2)
               << imt.msm(3) << imt.msm(4) << imt.msm(5) << imt.msm(6)
               << imt.msm(7) << imt.msm(8) << std::arg(imt.imt(2)) << NAN
               << std::endl;
        }
    } else if (mode == "scan_angle") {
        const double ellip = std::sqrt(1 - fsq(aspect));
        Datafile df(std::cout << "# angle area perimeter q2 q3 q4\n"
                              << "# aspect ratio " << aspect << " ellipticity "
                              << ellip << "\n");
        for (auto angle : linspace(0, .5 * TWO_PI, 50, false)) {
            (ds.*method)(Sersic(1., aspect, angle));
            if (!interpolated_marching_squares)
                imt = imt_regular_marching_squares(ds, threshold);
            else
                imt = imt_interpolated_marching_squares(ds, threshold);
            df << angle << imt.area() << imt.perimeter() << imt.msm(2)
               << imt.msm(3) << imt.msm(4) << imt.msm(5) << imt.msm(6)
               << imt.msm(7) << imt.msm(8) << std::arg(imt.imt(2)) << NAN
               << std::endl;
        }
    }

    return 0;
}
