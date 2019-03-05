// imganalysis: trivial program which load a PNG file, and computes
// IMT's for a number of interpolated marching squares thresholds
// using the preliminary Papaya2 library.
// this program considers only the red color channel.
// 2019 Sebastian Kapfer <sebastian.kapfer@fau.de>
#include "papaya2.hpp"
#include "picopng.hpp"
#include "readarg.hpp"

using namespace papaya2;

using bytevec_t = std::vector<unsigned char>;

static void load_file(bytevec_t *buffer, string const &filename);

struct PingFile
{
    unsigned long width_, height_;
    bytevec_t buffer_;

    enum { RED, GREEN, BLUE, ALPHA };

    PingFile(const string &infilename)
    {
        bytevec_t png_file_data;
        load_file(&png_file_data, infilename);
        if (!png_file_data.size())
            throw std::runtime_error("file is empty");
        int error = decodePNG(buffer_, width_, height_, png_file_data.data(),
                              png_file_data.size(), true);
        if (error != 0)
            throw std::runtime_error("picopng error " + std::to_string(error));
        if (buffer_.size() != width_ * height_ * 4)
            throw std::runtime_error("picopng didn't decode right");
    }

    int operator()(int i, int j) const
    {
        if (i == -1 || j == -1)
            return 0;
        if (i == width() || j == height())
            return 0;
        return const_cast<PingFile *>(this)->at(i, j, RED);
    }

    unsigned char &at(int i, int j, int c)
    {
        if (!(i >= 0 && i < width()))
            throw std::range_error("PingFile i");
        if (!(j >= 0 && j < height()))
            throw std::range_error("PingFile j");
        if (!(c >= 0 && c < 4))
            throw std::range_error("PingFile c");
        j = height() - 1 - j;
        int index = (i + j * width()) * 4 + c;
        return buffer_.at(index);
    }

    int width() const { return width_; }
    int height() const { return height_; }
    double pixel_width() const { return 1.; }
    double pixel_height() const { return 1.; }
};

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

void load_file(bytevec_t *buffer, const string &filename)
{
    std::ifstream file(filename.c_str(),
                       std::ios::in | std::ios::binary | std::ios::ate);

    // get filesize
    std::streamsize size = 0;
    if (file.seekg(0, std::ios::end).good())
        size = file.tellg();
    if (file.seekg(0, std::ios::beg).good())
        size -= file.tellg();

    // read contents of the file into the vector
    if (size > 0) {
        buffer->resize((size_t)size);
        file.read((char *)buffer->data(), size);
    } else
        buffer->clear();
}
