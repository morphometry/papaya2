#pragma once

// Papaya2
// header-only library for computing 2D irreducible Minkowski tensors
// 2018-2019 Sebastian Kapfer <sebastian.kapfer@fau.de>
// 2019 Fabian Schaller <physik@fabian-schaller.de>

#include "tools.hpp"
#include <limits>

namespace papaya2 {

// BasicPhoto
//
// This class is an example for the (static) interface for all
// PHOTO types that are intended to be used as image sources.
// should provide the proper member functions, though they
// need not necessarily derive from Photo.
//
// required functions for a PHOTO are:
//
// operator() for pixel access.
//
// width() and height() return the width and height, in pixels.
// please return int, not unsigned int.
//
// The image is considered to consist of the pixels with indices i,j
// where 0 <= i < width() and 0 <= height().
//
// You can use make_padded_view to supply padding pixels around an
// image which doesn't have any, s.t. any contours are closed.
//
// the functions origin() and upper_right() are used to convert to
// physical units.
template <typename TYPE> struct BasicPhoto
{
    void set_coordinates(double x0, double y0, double x1, double y1,
                         int width /* number of pixels */,
                         int height /* dito */)
    {
        origin_ = {x0, y0};
        ur_ = {x1, y1};
        width_ = width;
        height_ = height;
        lenx = (x1 - x0) / width_;
        leny = (y1 - y0) / height_;
        data.resize(size_t(width + 2) * size_t(height + 2));
    }

    // helper method to fill the Photo with a discretized
    // (sampled) version of the provided function
    // use set_coordinates first to define the coordinates.
    // the method below uses integration instead of sampling.
    template <typename FUNCTION> void sample_function(const FUNCTION &function)
    {
        for (int j = 0; j < height(); ++j) {
            for (int i = 0; i < width(); ++i) {
                double x = (i + .5) * lenx + origin_[0];
                double y = (j + .5) * leny + origin_[1];
                at(i, j) = function(x, y);
            }
        }
    }

    // trapezoidal rule to numerically integrate the function over
    // the pixel area.
    template <typename FUNCTION>
    void integrate_function(const FUNCTION &function)
    {
        for (int j = 0; j < height(); ++j)
            for (int i = 0; i < width(); ++i) {
                const int ss = 5;
                int total_w = 0;
                for (int subj = 0; subj <= ss; ++subj)
                    for (int subi = 0; subi <= ss; ++subi) {
                        double x = (i + subi * 1. / ss) * lenx + origin_[0];
                        double y = (j + subj * 1. / ss) * leny + origin_[1];
                        int w = (1 + (subi != 0 && subi != ss)) *
                                (1 + (subj != 0 && subj != ss));
                        total_w += w;
                        at(i, j) += w * function(x, y);
                    }
                at(i, j) /= total_w;
            }
    }

    // read access to the image
    TYPE operator()(int i, int j) const { return at(i, j); }

    // write access
    TYPE &operator()(int i, int j) { return at(i, j); }

    int width() const { return width_; }
    int height() const { return height_; }
    vec_t origin() const { return origin_; }
    vec_t upper_right() const { return ur_; }

  protected:
    int width_, height_; // number of pixels
    double lenx, leny;   // dimensions of each pixel
    vec_t origin_, ur_;
    std::vector<TYPE> data;

    size_t pixel_index(int i, int j) const
    {
        if (i < 0 || i >= width() || j < 0 || j >= height())
            throw std::range_error("invalid pixel index in BasicPhoto");
        return size_t(j + 1) * size_t(width() + 1) + size_t(i) + 1u;
    }

    TYPE &at(int i, int j) { return data[pixel_index(i, j)]; }

    const TYPE &at(int i, int j) const { return data[pixel_index(i, j)]; }
};

// compute pixel diagonal from coordinate system
template <typename PHOTO> vec_t pixel_diagonal(PHOTO const &photo)
{
    vec_t const diagonal = photo.upper_right() - photo.origin();
    return {diagonal[0] / photo.width(), diagonal[1] / photo.height()};
}

using Photo = BasicPhoto<double>;

// helper class to define transformed photos.
template <typename PHOTO> struct PhotoAdapter
{
    const PHOTO &original;

    PhotoAdapter(const PHOTO &ph) : original(ph) {}

    int width() const { return original.width(); }
    int height() const { return original.height(); }
    vec_t origin() const { return original.origin(); }
    vec_t upper_right() const { return original.upper_right(); }
};

// adapter to threshold a photo.
// returns photo >= threshold, i.e. ones and zeros.
// do not use directly, use make_thresholded_view.
template <typename PHOTO, typename THRESHOLD>
struct ThresholdingAdapter : PhotoAdapter<PHOTO>
{
    using PhotoAdapter<PHOTO>::original;
    const THRESHOLD threshold;

    ThresholdingAdapter(const PHOTO &ph, const THRESHOLD &t)
        : PhotoAdapter<PHOTO>(ph), threshold(t)
    {}

    int operator()(int x, int y) const
    {
        return int(original(x, y) >= threshold);
    }
};

template <typename PHOTO, typename THRESHOLD>
auto make_thresholded_view(const PHOTO &p, const THRESHOLD &t)
    -> ThresholdingAdapter<PHOTO, THRESHOLD>
{
    return ThresholdingAdapter<PHOTO, THRESHOLD>(p, t);
}

// adapter to pad a photo.
// returns photo(i,j) if (i,j) is in the original image,
// and padding_value otherwise.  data in the photo is cast to double.
// do not use directly, use make_padded_view.
template <typename PHOTO> struct PaddingAdapter : PhotoAdapter<PHOTO>
{
    using PhotoAdapter<PHOTO>::original;
    double const padding_value;

    PaddingAdapter(const PHOTO &ph, double pv)
        : PhotoAdapter<PHOTO>(ph), padding_value(pv)
    {}

    double operator()(int i, int j) const
    {
        if (i == 0 || j == 0)
            return padding_value;
        if (i == original.width() + 1 || j == original.height() + 1)
            return padding_value;
        return double(original(i - 1, j - 1));
    }

    int width() const { return original.width() + 2; }
    int height() const { return original.height() + 2; }
    vec_t origin() const
    {
        return original.origin() - pixel_diagonal(original);
    }
    vec_t upper_right() const
    {
        return original.upper_right() + pixel_diagonal(original);
    }
};

template <typename PHOTO>
vec_t pixel_diagonal(const PaddingAdapter<PHOTO> &photo)
{
    // padded photo has same-size pixels
    return pixel_diagonal(photo.original);
}

template <typename PHOTO>
auto make_padded_view(const PHOTO &p, double pv) -> PaddingAdapter<PHOTO>
{
    return PaddingAdapter<PHOTO>(p, pv);
}

template <typename TYPE, typename PHOTO>
void find_min_max(TYPE *min, TYPE *max, const PHOTO &photo)
{
    int const width = photo.width();
    int const height = photo.height();

    TYPE n = std::numeric_limits<TYPE>::max();
    TYPE x = std::numeric_limits<TYPE>::lowest();

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            n = std::min(n, TYPE(photo(i, j)));
            x = std::min(x, TYPE(photo(i, j)));
        }
    }

    if (min)
        *min = n;
    if (max)
        *max = x;
}

// result (return type) of IMT computations.
// holds the values of the IMT's and, in addition, the area.
// FIXME proper accumulators
struct MinkowskiAccumulator
{
    // for now, the upper limit in S is statically defined here, as the main
    // cost seems to be in the geometry processing, and a few extra cexp's
    // don't cost as much.  we leave it to the user to ignore any data she
    // does not need.
    static const int MAX_S = 12;

    MinkowskiAccumulator()
    {
        area_ = peri_ = 0.;
        for (int s : range(2, MAX_S + 1))
            psi(s) = 0.;
    }

    // add areas.
    friend void add_triangle_area(MinkowskiAccumulator *acc, vec_t, vec_t v0,
                                  vec_t v1, vec_t v2)
    {
        acc->area_ += kahan_triangle_area(v0, v1, v2);
    }

    template <typename CONTAINER>
    friend void add_polygon_area(MinkowskiAccumulator *acc,
                                 const CONTAINER &vertices)
    {
        acc->area_ += shoelace_formula(vertices);
    }

    // add a piece of straight-edge contour.
    friend void add_contour_segment(MinkowskiAccumulator *acc, vec_t offset,
                                    vec_t begin, vec_t end)
    {
        (void)offset;
        double const p = std::hypot(end[0] - begin[0], end[1] - begin[1]);
        if (p == 0.)
            return;
        complex_t n = {(end[1] - begin[1]) / p, -(end[0] - begin[0]) / p};
        acc->peri_ += p;
        for (int s : range(2, MAX_S + 1))
            acc->psi(s) += p * std::pow(n, s);
    }

    // accessors to read data
    double area() const { return area_; }
    double perimeter() const { return peri_; }

    double msm(int s) const
    {
        s = std::abs(s);
        if (s >= 2 && s <= MAX_S)
            return std::abs(psi(s)) / peri_;
        else
            throw std::logic_error("msm(" + std::to_string(s) + ") called");
    }

    complex_t imt(int s) const
    {
        if (s == 0)
            return perimeter();
        else if (s < 0)
            return std::conj(imt(-s));
        else if (s >= 2 && s <= MAX_S)
            return psi(s);
        else
            throw std::logic_error("imt(" + std::to_string(s) + ") called");
    }

    double beta102() const { return (1. - msm(2)) / (1. + msm(2)); }
    // isoperimetric ratio
    double isoper() const { return 2 * TWO_PI * area() / fsq(perimeter()); }

  private:
    static double kahan_triangle_area(vec_t v0, vec_t v1, vec_t v2)
    {
        double a = norm(v0 - v1);
        double b = norm(v1 - v2);
        double c = norm(v2 - v0);
        // sort
        if (a < b)
            std::swap(a, b);
        if (b < c)
            std::swap(b, c);
        if (a < b)
            std::swap(a, b);
        double x =
            (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c));
        return .25 * std::sqrt(x);
    }

    template <typename CONTAINER>
    static double shoelace_formula(const CONTAINER &vertices)
    {
        const int N = vertices.size();
        double twice_area = 0.;
        // shoelace formula for the area
        for (int i = 0; i != N; ++i) {
            auto const &begin = vertices.at(i);
            auto const &end = vertices.at((i + 1) % N);
            twice_area += begin[0] * end[1] - end[0] * begin[1];
        }
        return .5 * twice_area;
    }

    double area_, peri_;
    complex_t &psi(int s) { return int_psi_.at(s - 2); }
    complex_t psi(int s) const { return int_psi_.at(s - 2); }
    std::array<complex_t, MAX_S - 1> int_psi_;
};

// utility class which can be used with marching squares to
// write an isocontour to a file which may be plotted
// by plot "..." w vec
struct GnuplottableContour
{
    GnuplottableContour(const string &filename) : df(filename) {}
    GnuplottableContour(std::ostream &os) : df(os) {}

    friend void add_contour_segment(GnuplottableContour *sink, vec_t off,
                                    vec_t begin, vec_t end)
    {
        sink->df << (off[0] + begin[0]) << (off[1] + begin[1])
                 << (end[0] - begin[0]) << (end[1] - begin[1]) << std::endl;
    }

  private:
    Datafile df;
};

// add a polygonal contour.
// by default, decompose it into edges.
template <typename SINK, typename CONTAINER>
void add_polygon_contour(SINK *sink, const CONTAINER &vertices)
{
    const int N = vertices.size();
    for (int i = 0; i != N; ++i) {
        auto const &begin = vertices.at(i);
        auto const &end = vertices.at((i + 1) % N);
        add_contour_segment(sink, {0., 0.}, begin, end);
    }
}

// add a triangular area.
// by default, do nothing.
template <typename SINK>
void add_triangle_area(SINK *, vec_t offset, vec_t v0, vec_t v1, vec_t v2)
{
    (void)offset;
    (void)v0;
    (void)v1;
    (void)v2;
}

// add a polygonal area.
// by default, do nothing.
template <typename SINK, typename CONTAINER>
void add_polygon_area(SINK *sink, const CONTAINER &vertices)
{
    (void)sink;
    (void)vertices;
}

struct MarchingSquaresFlags
{
    explicit MarchingSquaresFlags(unsigned i = 0) : i_(i) {}
    using cref = const MarchingSquaresFlags &;

    friend bool operator&(cref lhs, cref rhs)
    {
        return (lhs.i_ & rhs.i_) != 0u;
    }

    MarchingSquaresFlags &operator|=(cref rhs)
    {
        i_ |= rhs.i_;
        return *this;
    }

    friend MarchingSquaresFlags operator|(cref lhs, cref rhs)
    {
        return MarchingSquaresFlags(lhs.i_ | rhs.i_);
    }

    MarchingSquaresFlags &operator^=(cref rhs)
    {
        i_ ^= rhs.i_;
        return *this;
    }

  private:
    unsigned i_;
};

static auto const ANALYZE_WHITE = MarchingSquaresFlags(0);
static auto const ANALYZE_BLACK = MarchingSquaresFlags(1);
static auto const CONNECT_WHITE = MarchingSquaresFlags(0);
static auto const CONNECT_BLACK = MarchingSquaresFlags(2);
// FIXME implement
// static auto const DISABLE_INTERPOLATION = MarchingSquaresFlags(4);

// interpolated marching squares, core routine handling a 2x2 neighborhood
// FIXME this is missing the curvature measures
template <typename SINK>
void add_interpolated_four_neighborhood(SINK *sink, vec_t const &off,
                                        vec_t const &pix_diag, double ll,
                                        double ul, double lr, double ur,
                                        double threshold,
                                        MarchingSquaresFlags flags)
{
    unsigned lut_index;

    if (flags & ANALYZE_BLACK)
    {
        lut_index = (ll < threshold) * 1 + (ul < threshold) * 2 +
                     (lr < threshold) * 4 + (ur < threshold) * 8;
        flags ^= CONNECT_BLACK;
    }
    else
    {
        lut_index = (ll >= threshold) * 1 + (ul >= threshold) * 2 +
                    (lr >= threshold) * 4 + (ur >= threshold) * 8;
    }

    // interpolate between 0.5 and 1.5.
    auto msq_interp = [](double left, double threshold,
                         double right) -> double {
        return (threshold - left) / (right - left) + .5;
    };

    // vertices, located on the edges of the neighborhood (NOT pixels)
    // i.e. on lines connecting pixel centers
    vec_t lower = {msq_interp(ll, threshold, lr), 0.5};
    vec_t upper = {msq_interp(ul, threshold, ur), 1.5};
    vec_t left = {0.5, msq_interp(ll, threshold, ul)};
    vec_t right = {1.5, msq_interp(lr, threshold, ur)};
    vec_t sw = {0.5, 0.5};
    vec_t se = {1.5, 0.5};
    vec_t nw = {0.5, 1.5};
    vec_t ne = {1.5, 1.5};
    // scale to physical coordinates
    lower = elementwise_product(lower, pix_diag);
    upper = elementwise_product(upper, pix_diag);
    left = elementwise_product(left, pix_diag);
    right = elementwise_product(right, pix_diag);
    sw = elementwise_product(sw, pix_diag);
    se = elementwise_product(se, pix_diag);
    nw = elementwise_product(nw, pix_diag);
    ne = elementwise_product(ne, pix_diag);
    // look-up table
    switch (lut_index) {
    case 0:
        // no area
        // no perimeter
        break;
    case 1:
        add_triangle_area(sink, off, lower, left, sw); // sw corner
        add_contour_segment(sink, off, lower, left);
        break;
    case 2:
        add_triangle_area(sink, off, left, upper, nw); // nw corner
        add_contour_segment(sink, off, left, upper);
        break;
    case 3:
        add_triangle_area(sink, off, lower, upper, nw); // left half
        add_triangle_area(sink, off, nw, sw, lower);
        add_contour_segment(sink, off, lower, upper);
        break;
    case 4:
        add_triangle_area(sink, off, right, lower, se); // se half
        add_contour_segment(sink, off, right, lower);
        break;
    case 5:
        add_triangle_area(sink, off, right, left, sw); // lower half
        add_triangle_area(sink, off, sw, se, right);
        add_contour_segment(sink, off, right, left);
        break;
    case 6:
        if (flags & CONNECT_BLACK) {
            add_contour_segment(sink, off, left, upper);
            add_triangle_area(sink, off, nw, left, upper);
            add_contour_segment(sink, off, right, lower);
            add_triangle_area(sink, off, lower, se, right);
        } else {
            add_triangle_area(sink, off, nw, left, upper);
            add_triangle_area(sink, off, left, lower, upper);
            add_triangle_area(sink, off, upper, lower, right);
            add_triangle_area(sink, off, lower, se, right);
            add_contour_segment(sink, off, left, lower);
            add_contour_segment(sink, off, right, upper);
        }
        break;
    case 7:
        add_triangle_area(sink, off, upper, nw, sw); // missing ne corner
        add_triangle_area(sink, off, right, upper, sw);
        add_triangle_area(sink, off, se, right, sw);
        add_contour_segment(sink, off, right, upper);
        break;
    case 8:
        add_triangle_area(sink, off, upper, right, ne); // ne corner
        add_contour_segment(sink, off, upper, right);
        break;
    case 9:
        if (flags & CONNECT_BLACK) {
            add_contour_segment(sink, off, lower, left);
            add_triangle_area(sink, off, sw, lower, left);
            add_contour_segment(sink, off, upper, right);
            add_triangle_area(sink, off, upper, right, ne);
        } else {
            add_triangle_area(sink, off, sw, lower, left);
            add_triangle_area(sink, off, lower, right, left);
            add_triangle_area(sink, off, left, right, upper);
            add_triangle_area(sink, off, upper, right, ne);
            add_contour_segment(sink, off, lower, right);
            add_contour_segment(sink, off, upper, left);
        }
        break;
    case 10:
        add_triangle_area(sink, off, nw, left, right); // upper half
        add_triangle_area(sink, off, right, ne, nw);
        add_contour_segment(sink, off, left, right);
        break;
    case 11:
        add_triangle_area(sink, off, right, ne, nw); // missing se corner
        add_triangle_area(sink, off, lower, right, nw);
        add_triangle_area(sink, off, sw, lower, nw);
        add_contour_segment(sink, off, lower, right);
        break;
    case 12:
        add_triangle_area(sink, off, upper, lower, se); // right half
        add_triangle_area(sink, off, se, ne, upper);
        add_contour_segment(sink, off, upper, lower);
        break;
    case 13:
        add_triangle_area(sink, off, ne, upper, se); // missing nw corner
        add_triangle_area(sink, off, upper, left, se);
        add_triangle_area(sink, off, left, sw, se);
        add_contour_segment(sink, off, upper, left);
        break;
    case 14:
        add_triangle_area(sink, off, nw, left, ne); // missing sw corner
        add_triangle_area(sink, off, left, lower, ne);
        add_triangle_area(sink, off, lower, se, ne);
        add_contour_segment(sink, off, left, lower);
        break;
    case 15:
        add_triangle_area(sink, off, sw, se, ne); // full square
        add_triangle_area(sink, off, ne, nw, sw);
        // no perimeter
        break;
    default:
        std::abort();
    }
}

// interpolated marching squares, loop over the whole image
// effectively calls add_area and add_perimeter on "sink" once
// a piece of the contour has been identified.
// processes the whole photo, (w-2)x(h-2) neighborhoods,
template <typename SINK, typename PHOTO, typename THRESHOLD>
void trace_isocontour_interpolated_marching_squares(
    SINK *sink, const PHOTO &ph, const THRESHOLD &threshold,
    MarchingSquaresFlags flags = MarchingSquaresFlags())
{
    const int width = ph.width();
    const int height = ph.height();

    vec_t const pix_diag = pixel_diagonal(ph);
    vec_t const origin = ph.origin();

    for (int j = 0; j < height - 1; ++j)
        for (int i = 0; i < width - 1; ++i) {
            vec_t off = origin + vec_t{i * pix_diag[0], j * pix_diag[1]};
            add_interpolated_four_neighborhood(
                sink, off, pix_diag, ph(i, j), ph(i, j + 1), ph(i + 1, j),
                ph(i + 1, j + 1), threshold, flags);
        }
}

// convenience wrapper to compute IMT's with interpolated marching squares
template <typename PHOTO>
MinkowskiAccumulator imt_interpolated_marching_squares(
    const PHOTO &ph, double threshold,
    MarchingSquaresFlags flags = MarchingSquaresFlags())
{
    MinkowskiAccumulator acc;
    trace_isocontour_interpolated_marching_squares(&acc, ph, threshold, flags);
    return acc;
}

// convenience wrapper to compute IMT's with regular marching squares
// reproduce (inferior) results of regular marching squares
// by applying interpolated marching squares on a binarized image
template <typename PHOTO>
MinkowskiAccumulator imt_regular_marching_squares(
    const PHOTO &ph, double threshold,
    MarchingSquaresFlags flags = MarchingSquaresFlags())
{
    MinkowskiAccumulator acc;
    auto thr_ph = make_thresholded_view(ph, threshold);
    trace_isocontour_interpolated_marching_squares(&acc, thr_ph, .5, flags);
    return acc;
}

// compute the IMTs of a single polygon
template <typename CONTAINER>
MinkowskiAccumulator imt_polygon(const CONTAINER &vertices)
{
    MinkowskiAccumulator acc;
    add_polygon_area(&acc, vertices);
    add_polygon_contour(&acc, vertices);
    return acc;
}

using complex_image_t = BasicPhoto<complex_t>;

// FIXME expose flags
template <typename PHOTO, typename THRESHOLD>
void minkowski_map_interpolated_marching_squares(complex_image_t *out,
                                                 const PHOTO &ph,
                                                 const THRESHOLD &threshold,
                                                 int s)
{
    const int lastx = ph.width() - 2;
    const int lasty = ph.height() - 2;

    vec_t const origin = ph.origin();
    vec_t const pix_diag = pixel_diagonal(ph);
    vec_t const half_a_pixdiag = pix_diag / 2;
    vec_t mmap_origin = ph.origin() + half_a_pixdiag;
    vec_t mmap_upperright = ph.upper_right() - half_a_pixdiag;
    out->set_coordinates(mmap_origin[0], mmap_origin[1], mmap_upperright[0],
                         mmap_upperright[1], lastx + 1, lasty + 1);

    for (int j = 0; j <= lasty; ++j) {
        for (int i = 0; i <= lastx; ++i) {
            vec_t const off = origin + vec_t{i * pix_diag[0], j * pix_diag[1]};
            MinkowskiAccumulator minkval;
            add_interpolated_four_neighborhood(
                &minkval, off, pix_diag, ph(i, j), ph(i, j + 1), ph(i + 1, j),
                ph(i + 1, j + 1), threshold, MarchingSquaresFlags());
            (*out)(i, j) = minkval.imt(s);
        }
    }
}

} // namespace papaya2
