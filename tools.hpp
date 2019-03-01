#pragma once

#include "vector.hpp"
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace papaya2 {

using complex_t = std::complex<double>;
using string = std::string;
using vec_t = std::array<double, 2>;
using namespace vector_math_for_std_array;

static double const TWO_PI = 6.2831853071795864769252867665590057683943;
static vec_t const ORIGIN = {0., 0.};
static double const SQRT2 = std::sqrt(2);

inline double atan2(const vec_t &v) { return std::atan2(v[1], v[0]); }

inline vec_t cos_sin(double phi) { return {std::cos(phi), std::sin(phi)}; }

inline double fsq(double x) { return x * x; }

inline double deg2rad(double phi) { return phi * (TWO_PI / 360.); }

// produce scales with evenly and log-evenly spaced datapoints.
// set last argument to true if you want "high" to be the last datapoint.
inline std::vector<double> logspace(double low, double high, int n = 10,
                                    bool include_endpoint = false)
{
    double logincr = (log(high) - log(low)) / (n - include_endpoint);
    std::vector<double> ret;
    ret.reserve(n);
    for (int i = 0; i != n; ++i)
        ret.push_back(exp(i * logincr + log(low)));
    return ret;
}

inline std::vector<double> linspace(double low, double high, int n = 10,
                                    bool include_endpoint = false)
{
    double incr = (high - low) / (n - include_endpoint);
    std::vector<double> ret;
    ret.reserve(n);
    for (int i = 0; i != n; ++i)
        ret.push_back(i * incr + low);
    return ret;
}

// container which contains a sequence of integers
struct IntegerRange
{
    using value_type = int;

    struct iterator
    {
        using difference_type = int;
        using value_type = int;
        using pointer = void;
        using reference = void;
        using iterator_category = std::bidirectional_iterator_tag;

        iterator() : x_(42) {}
        iterator(int x) : x_(x) {}

        int operator*() const { return x_; }

        iterator &operator++()
        {
            ++x_;
            return *this;
        }

        iterator &operator--()
        {
            --x_;
            return *this;
        }

        friend bool operator!=(const iterator &lhs, const iterator &rhs)
        {
            return lhs.x_ != rhs.x_;
        }

      private:
        int x_;
    };

    iterator begin() const { return min_; }
    iterator end() const { return max_excluded_; }

    int min_, max_excluded_;
};

inline IntegerRange range(int min, int max_excluded)
{
    return {min, max_excluded};
}
inline IntegerRange range(/* min*/ int max_excluded)
{
    return {0, max_excluded};
}

// auxiliary convenience class to write space-separated datafiles.
// for an example use see banana.cpp.
struct Datafile
{
    using ostream = std::ostream;
    using iomanip_t = ostream &(*)(ostream &);

    Datafile(ostream &os_)
    {
        // re-use the fstream as a regular stream, with a
        // foreign streambuf
        static_cast<ostream &>(os).rdbuf(os_.rdbuf());
        init();
    }

    Datafile(const string &filename) : os(filename) { init(); }

    friend Datafile &operator<<(Datafile &df, double x)
    {
        df.os << x << ' ';
        return df;
    }

    // hijack the std::endl manipulator to mean `end of data line'.
    // we don't actually flush, merely insert a newline.
    friend Datafile &operator<<(Datafile &df, iomanip_t x)
    {
        if (static_cast<iomanip_t>(&std::endl) == x) {
            // insert a newline.  (don't flush.)
            df.os << '\n';
            return df;
        } else if (static_cast<iomanip_t>(&std::flush) == x) {
            // flush
            df.os.flush();
            return df;
        }
        throw std::runtime_error("illegal manipulator inserted into Datafile");
    }

    void comment(const string &comment) { os << "# " << comment << "\n"; }

  private:
    std::ofstream os;
    void init()
    {
        os << std::scientific;
        os.precision(15);
    }
};

} // namespace papaya2
