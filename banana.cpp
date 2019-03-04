// banana: trivial program which load a FITS file, and computes
// IMT's for a number of interpolated marching squares thresholds
// using the preliminary Papaya2 library.
// 2018 Sebastian Kapfer <sebastian.kapfer@fau.de>
// 2018 Jenny Wagner <j.wagner@uni-heidelberg.de>
#include "papaya2.hpp"
#include "readarg.hpp"
#ifndef HAVE_CCFITS
#error banana requires FITS support to build (see Makefile)
#endif
#include <CCfits/CCfits>
#include <valarray>

struct FitsFile : Photo
{
    FitsFile(const string &infilename, unsigned hdu = 0)
    {
        if (hdu != 0)
            throw std::runtime_error("hdu != 0 not implemented");
        CCfits::FITS::setVerboseMode(true);
        CCfits::FITS infile(infilename.c_str(), CCfits::Read, true);
        CCfits::PHDU &phdu = infile.pHDU();
        // FIXME extract metadata from FITS to convert pixel coordinates
        // to sky coordinates.  for now, we map the entire image
        // to the unit square.
        // FIXME readAllKeys often throws a warning of the type
        //       Fits Error: Keyword not found: HISTORY                                                                                                  
        //       Fits Error: Keyword not found: COMMENT                                                                                                  
        // at this point, which is harmless.  find a way to disable it.
        phdu.readAllKeys();
        set_coordinates(0, 0, 1, 1, phdu.axis(0), phdu.axis(1));
        std::valarray<double> buffer;
        phdu.read(buffer);
        const int w = width();
        const int h = height();
        double max = -1e99;
        double min = 1e99;
        for (int i = 0; i < w; ++i)
            for (int j = 0; j < h; ++j) {
                at(i, j) = buffer[(h - j - 1) * w + i];
                max = fmax(max, at(i, j));
                min = fmin(min, at(i, j));
            }

        std::cerr << "info: loaded " << w << "x" << h
                  << " FITS file, max lumin = " << max << ", min = " << min
                  << "\n";
    }

    // allow overwriting pixels by the user
    double &at(int i, int j) { return Photo::at(i, j); }

    void addnoise() { at(2, 2) = 1; }
};

int main(int argc, const char **argv)
{
    string infilename, outfilename;
    string maskfilename;
    int num_thresh = 10;
    double min_thresh = .01;
    double max_thresh = 1.;

    // process command-line arguments
    for (++argv; *argv; ++argv) {
        if (string(*argv) == "in") {
            infilename = read_arg<string>(argv++);
        } else if (string(*argv) == "mask") {
            maskfilename = read_arg<string>(argv++);
        } else if (string(*argv) == "out") {
            outfilename = read_arg<string>(argv++);
        } else if (string(*argv) == "mint") {
            min_thresh = read_arg<double>(argv++);
        } else if (string(*argv) == "maxt") {
            max_thresh = read_arg<double>(argv++);
        } else if (string(*argv) == "numt") {
            num_thresh = read_arg<unsigned long>(argv++);
        } else {
            std::cerr << "illegal argument: " << *argv << "\n";
            return 1;
        }
    }

    FitsFile infile(infilename);

    if (maskfilename != "") {
        // apply a mask to the image to extract a particular object in the
        // image.
        FitsFile mask(maskfilename);
        const int width = infile.width();
        const int height = infile.height();
        if (mask.width() != width || mask.height() != height) {
            std::cerr << "dimensions of image and mask don't match: " << width
                      << "x" << height << " vs. " << mask.width() << "x"
                      << mask.height() << ". aborting.\n";
            return 1;
        }

        int mini = width, maxi = -1;
        int minj = height, maxj = -1;

        for (int j = 0; j < height - 1; ++j)
            for (int i = 0; i < width - 1; ++i) {
                const double w = mask.at(i, j);
                if (!(w >= 0 && w <= 1.)) {
                    std::cerr << "invalid mask value: " << w << " at " << i
                              << ", " << j << "\n";
                    return 1;
                }
                if (w != 0.) {
                    mini = std::min(i, mini);
                    maxi = std::max(i, maxi);
                    minj = std::min(j, minj);
                    maxj = std::max(j, maxj);
                }
                infile.at(i, j) *= w;
            }

        std::cerr << "object bounding box lower-left corner (inclusive): "
                  << mini << " " << minj << "\n";
        std::cerr << "object bounding box upper-right corner (exclusive): "
                  << (maxi + 1) << " " << (maxj + 1) << "\n";
    }

#if 0
    // set a single pixel to 1
    infile.addnoise ();
#endif

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

#if 0
    // dump a contour for demonstration
    // gnuplot it via
    //          plot "contour.dat" w vec
    GnuplottableContour gc ("contour.dat");
    trace_isocontour_interpolated_marching_squares (&gc, infile, 1e-3);
#endif

    return 0;
}
