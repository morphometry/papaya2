// imganalysis: trivial program which load a PNG file, and computes
// IMT's for a number of interpolated marching squares thresholds
// using the preliminary Papaya2 library.
// this program considers only the red color channel.
// 2019 Sebastian Kapfer <sebastian.kapfer@fau.de>
#include <papaya2.hpp>
#include "readarg.hpp"

// pull in loading of PNG image files.
#include <papaya2/pingfile.hpp>

using namespace papaya2;

int main(int argc, const char **argv)
{
    string infilename, outfilename, contours_filename;
    int num_thresh = 10;
    double min_thresh = 1.;
    double max_thresh = 255.;

    // process command-line arguments
    for (++argv; *argv; ++argv) {
        if (string(*argv) == "in") {
            infilename = read_arg<string>(argv++);
        } else if (string(*argv) == "out") {
            outfilename = read_arg<string>(argv++);
        } else if (string(*argv) == "mint") {
            min_thresh = read_arg<double>(argv++);
        } else if (string(*argv) == "maxt") {
            max_thresh = read_arg<double>(argv++);
        } else if (string(*argv) == "numt") {
            num_thresh = read_arg<unsigned long>(argv++);
        } else if (string(*argv) == "contours") {
            contours_filename = read_arg<string>(argv++);
        } else if (string(*argv) == "help" || string(*argv) == "-h") {
            std::cerr << "for help, please see "
                         "https://morphometry.org/software/papaya2/\n";
            return 0;
        } else {
            std::cerr << "illegal argument: " << *argv << "\n";
            return 1;
        }
    }

    if (outfilename == "") {
        std::cerr << "output filename not set (use the 'out' argument)\n";
        return 1;
    }

    PingFile infile(infilename);
    Datafile outfile(outfilename);
    outfile.comment("threshold area perim q2 arg2 q3 arg3 q4 arg4 q5 arg5 q6 "
                    "arg6 q7 arg7 q8 arg8");

    for (auto thresh : logspace(min_thresh, max_thresh, num_thresh, true)) {
        auto imt = imt_interpolated_marching_squares(infile, thresh);
        outfile << thresh << imt.area() << imt.perimeter() << imt.msm(2)
                << std::arg(imt.imt(2)) << imt.msm(3) << std::arg(imt.imt(3))
                << imt.msm(4) << std::arg(imt.imt(4)) << imt.msm(5)
                << std::arg(imt.imt(5)) << imt.msm(6) << std::arg(imt.imt(6))
                << imt.msm(7) << std::arg(imt.imt(7)) << imt.msm(8)
                << std::arg(imt.imt(8)) << std::endl;
    }

    // dump a contour for demonstration
    // gnuplot it via
    //          plot "contour.dat" w vec
    if (contours_filename != "") {
        std::ofstream contours_file(contours_filename);
        GnuplottableContour gc(contours_file);
        for (auto thresh : logspace(min_thresh, max_thresh, num_thresh, true)) {
            trace_isocontour_interpolated_marching_squares(&gc, infile, thresh);
            contours_file << "\n";
        }
    }

    return 0;
}
