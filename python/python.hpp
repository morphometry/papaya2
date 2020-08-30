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

struct PythonErrorSet
{};

struct UniquePyPtr : std::unique_ptr<PyObject, RefDeccer>
{
    using base_t = std::unique_ptr<PyObject, RefDeccer>;
    using base_t::base_t;

    template <typename TYPE> TYPE *reinterpret() const
    {
        return reinterpret_cast<TYPE *>(get());
    }
};

} // namespace python
} // namespace papaya2
