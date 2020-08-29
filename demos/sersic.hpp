#pragma once

#include "papaya2.hpp"

// Sérsic profile with the given parameters.
// may be fed into "sample_function" and "integrate_function"
// to obtain test images to analyze.
namespace papaya2 {

struct Sersic
{
    double a, b;
    int exponent;
    double cos, sin;

    Sersic(double a_ = 1., double b_ = 1., double angle_ = 0.,
           int exponent_ = 1)
        : a(a_), b(b_), exponent(exponent_), cos(std::cos(angle_)),
          sin(std::sin(angle_))
    {
        if (exponent > 5 || exponent < 1)
            throw std::runtime_error("unsupported Sérsic exponent");
    }

    double operator()(double x, double y) const
    {
        // transform to diagonal system
        double rot_x = cos * x + sin * y;
        double rot_y = cos * y - sin * x;
        double rsq = fsq(rot_x / a) + fsq(rot_y / b);
        double z = std::pow(rsq, 1. / (2 * exponent));
        return std::exp(-bn(exponent) * (z - 1.));
    }

  private:
    static double bn(int exponent)
    {
        double bn[] = {-1, 1, 3.67206075, 5.67016119, 7.66924944, 9.66871461};
        return bn[exponent];
    }
};

} // namespace papaya2
