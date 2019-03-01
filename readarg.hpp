#pragma once

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>

namespace papaya2 {

using string = std::string;

// processing of command-line arguments
// handles error signalling (via exceptions) and type conversion.
template <typename TYPE> TYPE read_arg(const char **argv);

static std::pair<string, string> int_pop_arg(const char **argv)
{
    assert(*argv);
    string kw = *argv++;
    if (!*argv)
        throw std::runtime_error("missing argument for kw " + kw);
    return std::make_pair(kw, string(*argv));
}

template <> inline string read_arg(const char **argv)
{
    string kw, value;
    std::tie(kw, value) = int_pop_arg(argv);
    return value;
}

template <> inline unsigned long read_arg(const char **argv)
{
    string kw, value;
    std::tie(kw, value) = int_pop_arg(argv);
    try {
        return std::stoul(value);
    } catch (...) {
        throw std::runtime_error("cannot convert value " + value +
                                 " for keyword " + kw);
    }
}

template <> inline double read_arg(const char **argv)
{
    string kw, value;
    std::tie(kw, value) = int_pop_arg(argv);
    try {
        return std::stod(value);
    } catch (...) {
        throw std::runtime_error("cannot convert value " + value +
                                 " for keyword " + kw);
    }
}

// print error message and exit
inline void die(const string &message)
{
    std::cerr << message << std::endl;
    std::exit(1);
}

} // namespace papaya2
