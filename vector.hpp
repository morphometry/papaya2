#ifndef POSTLHC_VECTOR_HPP_INCLUDED 
#define POSTLHC_VECTOR_HPP_INCLUDED 
// (c) 2018-2019 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
// a few small routines for doing fixed-length vector operations with std::array
// (extend as necessary)
// these routines are packed in their own namespaces so they can be enabled as
// required.

#include <array>
#include <cassert>
#include <cmath>
#include <string>
#include <type_traits>

namespace vector_math_for_std_array {
namespace detail {
// effectively duplicate std::size to avoid depending on C++17
template <typename VEC> size_t vector_dimension(const VEC &v)
{
    return v.size();
}

template <typename TYPE, size_t DIM>
constexpr size_t vector_dimension(const TYPE (&)[DIM])
{
    return DIM;
}

template <typename TYPE, typename TYPE2, size_t DIM>
using return_array =
    std::array<typename std::common_type<TYPE, TYPE2>::type, DIM>;

template <typename VEC> struct SequencePrinter
{
    SequencePrinter(const VEC &v_, const char *sep) : v(v_), separator(sep) {}

    inline friend std::ostream &operator<<(std::ostream &os,
                                           const SequencePrinter<VEC> &printer)
    {
        printer.print_to(os);
        return os;
    }

  private:
    void print_to(std::ostream &os) const
    {
        const size_t size = detail::vector_dimension(v);
        if (size > 0) {
            os << v[0];
            for (size_t i = 1; i != size; ++i)
                os << separator << v[i];
        }
    }

    const VEC &v;
    const std::string separator;
};
} // namespace detail

template <typename VEC>
auto comma_separated(const VEC &v) -> detail::SequencePrinter<VEC>
{
    return {v, ", "};
}

template <typename VEC>
auto space_separated(const VEC &v) -> detail::SequencePrinter<VEC>
{
    return {v, " "};
}

template <typename TYPE, size_t DIM>
std::array<TYPE, DIM> operator+(const std::array<TYPE, DIM> &lhs,
                                const std::array<TYPE, DIM> &rhs)
{
    std::array<TYPE, DIM> ret;
    for (size_t n = 0; n != DIM; ++n)
        ret[n] = lhs[n] + rhs[n];
    return ret;
}

template <typename TYPE, size_t DIM, typename VEC2>
std::array<TYPE, DIM> &operator+=(std::array<TYPE, DIM> &lhs, const VEC2 &rhs)
{
    for (size_t n = 0; n != DIM; ++n)
        lhs[n] += rhs[n];
    return lhs;
}

template <typename TYPE, size_t DIM>
std::array<TYPE, DIM> operator-(const std::array<TYPE, DIM> &lhs,
                                const std::array<TYPE, DIM> &rhs)
{
    std::array<TYPE, DIM> ret;
    for (size_t n = 0; n != DIM; ++n)
        ret[n] = lhs[n] - rhs[n];
    return ret;
}

template <typename TYPE, size_t DIM, typename VEC2>
std::array<TYPE, DIM> &operator-=(std::array<TYPE, DIM> &lhs, const VEC2 &rhs)
{
    for (size_t n = 0; n != DIM; ++n)
        lhs[n] -= rhs[n];
    return lhs;
}

// scalar multiply
template <typename TYPE, size_t DIM, typename TYPE2>
std::array<TYPE, DIM> &operator*=(std::array<TYPE, DIM> &lhs, const TYPE2 &rhs)
{
    for (size_t n = 0; n != DIM; ++n)
        lhs[n] *= rhs;
    return lhs;
}

template <typename TYPE, size_t DIM, typename TYPE2>
auto operator*(const std::array<TYPE, DIM> &lhs, const TYPE2 &rhs)
    -> detail::return_array<TYPE, TYPE2, DIM>
{
    detail::return_array<TYPE, TYPE2, DIM> ret;
    for (size_t n = 0; n != DIM; ++n)
        ret[n] = lhs[n] * rhs;
    return ret;
}

template <typename TYPE, size_t DIM, typename TYPE2>
auto operator*(const TYPE2 &lhs, const std::array<TYPE, DIM> &rhs)
    -> decltype(rhs * lhs)
{
    return rhs * lhs;
}

template <typename TYPE, size_t DIM, typename TYPE2>
auto operator/(const std::array<TYPE, DIM> &lhs, const TYPE2 &rhs)
    -> detail::return_array<TYPE, TYPE2, DIM>
{
    using ret_t = detail::return_array<TYPE, TYPE2, DIM>;
    assert(rhs != 0);
    const typename ret_t::value_type scale = 1. / rhs;
    ret_t ret;
    for (size_t n = 0; n != DIM; ++n)
        ret[n] = lhs[n] * scale;
    return ret;
}

template <typename VEC1, typename VEC2>
auto inner_product(const VEC1 &lhs, const VEC2 &rhs)
    -> decltype(lhs[0] * rhs[0])
{
    using type = decltype(lhs[0] * rhs[0]);
    using detail::vector_dimension;
    const size_t size = vector_dimension(lhs);
    assert(size == vector_dimension(rhs));
    if (size == 0)
        return type();
    type ret = lhs[0] * rhs[0];
    for (size_t i = 1; i != size; ++i)
        ret += lhs[i] * rhs[i];
    return ret;
}

template <typename TYPE, size_t DIM, typename VEC2>
auto elementwise_product(const std::array<TYPE, DIM> &lhs, const VEC2 &rhs)
    -> std::array<decltype(lhs[0] * rhs[0]), DIM>
{
    std::array<decltype(lhs[0] * rhs[0]), DIM> ret;
    for (size_t n = 0; n != DIM; ++n)
        ret[n] = lhs[n] * rhs[n];
    return ret;
}

template <typename TYPE, size_t DIM>
auto norm_sq(const std::array<TYPE, DIM> &v) -> decltype(inner_product(v, v))
{
    return inner_product(v, v);
}

template <size_t DIM> double norm(const std::array<double, DIM> &v)
{
    return std::sqrt(norm_sq(v));
}

template <typename TYPE, size_t DIM, typename VEC2>
bool operator==(const std::array<TYPE, DIM> &lhs, const VEC2 &rhs)
{
    using detail::vector_dimension;
    assert(DIM == vector_dimension(rhs));
    for (size_t n = 0; n != DIM; ++n) {
        if (!(lhs[n] == rhs[n]))
            return false;
    }

    return true;
}
} // namespace vector_math_for_std_array

namespace vector_printing_for_std_array {

template <typename TYPE, size_t DIM>
std::ostream &operator<<(std::ostream &os, const std::array<TYPE, DIM> &rhs)
{
    for (size_t n = 0; n != DIM - 1; ++n)
        os << rhs[n] << ' ';
    return os << rhs[DIM - 1];
}

} // namespace vector_printing_for_std_array

#endif /* POSTLHC_VECTOR_HPP_INCLUDED */
