#pragma once

namespace papaya2 {

struct CheckPositivelyOriented
{
    static bool ordered(double x, double y, double z)
    {
        return x <= y && y <= z;
    }

    void add_triangle_area(vec_t, vec_t v0, vec_t v1, vec_t v2)
    {
        vec_t b = (v0 + v1 + v2) / 3;
        const double a0 = atan2(v0 - b);
        const double a1 = atan2(v1 - b);
        const double a2 = atan2(v2 - b);
        if (!ordered(a0, a1, a2) && !ordered(a1, a2, a0) &&
            !ordered(a2, a0, a1))
            throw std::logic_error("triangle is not positively oriented");
    }

    friend void add_contour_segment(CheckPositivelyOriented *, vec_t, vec_t,
                                    vec_t)
    {
        // do nothing
    }
};

template <typename PHOTO, typename FLAGS>
void check_msq_positively_oriented(const PHOTO &ph, const FLAGS &flags)
{
    CheckPositivelyOriented chk;
    trace_isocontour_interpolated_marching_squares(&chk, ph, .5, flags);
}

struct CApprox
{
    CApprox(double real, double imag, double scale, double tolerance = 1e-12)
        : ref_value_(real, imag), scale_(scale), tolerance_(tolerance)
    {}

    // FIXME allow overriding scale

    CApprox tolerance(double tolerance)
    {
        return CApprox(ref_value_.real(), ref_value_.imag(), scale_, tolerance);
    }

    friend bool operator==(const complex_t &is, const CApprox &approx)
    {
        const double diff = std::abs(is - approx.ref_value_);
        return diff / approx.scale_ < approx.tolerance_;
    }

    friend CApprox operator*(double factor, const CApprox &approx)
    {
        return CApprox(approx.ref_value_.real() * factor,
                       approx.ref_value_.imag() * factor,
                       approx.scale_ * factor, approx.tolerance_);
    }

    friend CApprox operator*(const CApprox &approx, double factor)
    {
        return factor * approx;
    }

    friend std::ostream &operator<<(std::ostream &os, const CApprox &approx)
    {
        return os << approx.ref_value_;
    }

  private:
    complex_t ref_value_;
    double scale_;
    const double tolerance_;
};

inline CApprox approx(double real, double imag, double scale = 1.)
{
    return CApprox(real, imag, scale);
}

inline CApprox approx(const complex_t &c, double extra_im = 0.,
                      double scale = 1.)
{
    return CApprox(c.real(), c.imag() + extra_im, scale);
}

} // namespace papaya2
