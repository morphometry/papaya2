#pragma once

#include <Python.h>
#include <memory>
#include <string>

namespace papaya2 {
namespace python {

struct RefDeccer
{
    void operator()(PyObject *ref) { Py_XDECREF(ref); }
};

struct UniquePyPtr : std::unique_ptr<PyObject, RefDeccer>
{
    using base_t = std::unique_ptr<PyObject, RefDeccer>;
    using base_t::base_t;

    template <typename TYPE> TYPE *reinterpret() const
    {
        return reinterpret_cast<TYPE *>(get());
    }
};


inline string to_utf8_string(PyObject *unicode)
{
#ifdef PYTHON_3
    if(!PyUnicode_Check(unicode))
        throw std::logic_error("non-PyUnicode passed into to_utf8_string");
    Py_ssize_t size;
    char const *buf = PyUnicode_AsUTF8AndSize(unicode, &size);
    return string(buf, buf + size);
#else
    if(!PyBytes_Check(unicode))
        throw std::logic_error("non-PyBytes passed into to_utf8_string");
    return PyBytes_AS_STRING(unicode);
#endif
}

struct KwargsIterator
{
    using value_type = std::pair<string, PyObject *>;

    KwargsIterator &operator++()
    {
        PyObject *key = 0;
        if(PyDict_Next(dict_, &pos_, &key, &key_and_value_.second)) {
            key_and_value_.first = to_utf8_string(key);
        } else {
            pos_ = 0; // end
        }
        return *this;
    }

    value_type const &operator*() const
    {
        return key_and_value_;
    }

    friend
    bool operator==(KwargsIterator const &lhs, KwargsIterator const &rhs)
    {
        return lhs.pos_ == rhs.pos_;
    }

    friend
    bool operator!=(KwargsIterator const &lhs, KwargsIterator const &rhs)
    {
        return lhs.pos_ != rhs.pos_;
    }

    KwargsIterator(PyObject *dict)
        : dict_(dict), pos_(0), key_and_value_("", nullptr)
    {
    }

    PyObject *dict_;
    Py_ssize_t pos_;
    value_type key_and_value_;
private:
    KwargsIterator();
};

struct Kwargs
{
    Kwargs(PyObject *dict)
        : dict_(dict)
    {
    }

    KwargsIterator begin() const
    {
        auto iter = KwargsIterator(dict_);
        ++iter;
        return iter;
    }

    KwargsIterator end() const
    {
        return KwargsIterator(dict_);
    }

    PyObject *dict_;
};

} // namespace python
} // namespace papaya2
