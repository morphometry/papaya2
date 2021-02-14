// Python/Numpy interface.
// we're not using Boost::Python here as we export very few
// functions and Boost::Python::Numpy is not in Debian.
// 2019-2020 Sebastian Kapfer <sebastian.kapfer@fau.de>
// TODO: Extract common code between the two functions
// TODO: Add test coverage for the pointpattern analysis

#include <papaya2.hpp>
#ifdef HAVE_CGAL
// without CGAL, there is no Voronoi.
#include <papaya2/voronoi.hpp>
#endif
#define PY_SSIZE_T_CLEAN
#include "python.hpp"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

using namespace papaya2;
using namespace papaya2::python;

static void np_set_double(UniquePyPtr &obj, int i, double value)
{
    auto arr = obj.reinterpret<PyArrayObject>();
    *reinterpret_cast<double *>(PyArray_GETPTR1(arr, i)) = value;
}

static void np_set_complex(UniquePyPtr &obj, int i, complex_t value)
{
    auto arr = obj.reinterpret<PyArrayObject>();
    *reinterpret_cast<complex_t *>(PyArray_GETPTR1(arr, i)) = value;
}

#ifdef HAVE_CGAL
static double np_get_double(UniquePyPtr &obj, int i)
{
    auto arr = obj.reinterpret<PyArrayObject>();
    return *reinterpret_cast<double const *>(PyArray_GETPTR1(arr, i));
}
#endif // HAVE_CGAL

static double np_get_double(UniquePyPtr &obj, int i, int j)
{
    auto arr = obj.reinterpret<PyArrayObject>();
    return *reinterpret_cast<double const *>(PyArray_GETPTR2(arr, i, j));
}

// make a new numpy array, 1D with length N, filled with NaN's
static UniquePyPtr np_make_new_vector(int N, int datatype)
{
    npy_intp dims[1] = {N};
    auto ret = UniquePyPtr(PyArray_SimpleNew(1, dims, datatype));
    if (!ret)
        throw std::runtime_error("np_make_new_vector (1)");
    auto nan = UniquePyPtr(Py_BuildValue("d", NAN));
    int ec =
        PyArray_FillWithScalar(ret.reinterpret<PyArrayObject>(), nan.get());
    if (ec != 0)
        throw std::runtime_error("np_make_new_vector (2)");
    return ret;
}

namespace {
    struct MinkValReturnData
    {
        MinkValReturnData(int N, int MAX_S)
        {
            ref_area_data = np_make_new_vector(N, NPY_DOUBLE);
            ref_peri_data = np_make_new_vector(N, NPY_DOUBLE);
            ref_msm_data = std::vector<UniquePyPtr>(MAX_S + 1);
            ref_imt_data = std::vector<UniquePyPtr>(MAX_S + 1);
            for (int s = 2; s <= MAX_S; ++s) {
                ref_msm_data[s] = np_make_new_vector(N, NPY_DOUBLE);
                ref_imt_data[s] = np_make_new_vector(N, NPY_CDOUBLE);
            }
        }

        UniquePyPtr ref_area_data;
        UniquePyPtr ref_peri_data;
        std::vector<UniquePyPtr> ref_msm_data;
        std::vector<UniquePyPtr> ref_imt_data;

        void assign(int index, MinkowskiAccumulator const &);
        UniquePyPtr move_to_dict();
    };

    void MinkValReturnData::assign(int index, MinkowskiAccumulator const &minkval)
    {
        np_set_double(ref_area_data, index, minkval.area());
        np_set_double(ref_peri_data, index, minkval.perimeter());
        for (unsigned s = 2; s < ref_msm_data.size(); ++s) {
            np_set_double(ref_msm_data[s], index, minkval.msm(s));
            np_set_complex(ref_imt_data[s], index, minkval.imt(s));
        }
    }

    UniquePyPtr MinkValReturnData::move_to_dict() {
        auto ref_return_dict = UniquePyPtr(PyDict_New());
        auto move_item_to_dict = [&](string const &name, UniquePyPtr &ref) {
            (void)PyDict_SetItemString(ref_return_dict.get(), name.c_str(),
                                       ref.release());
        };
        move_item_to_dict("area", ref_area_data);
        move_item_to_dict("perimeter", ref_peri_data);
        for (unsigned s = 2; s < ref_msm_data.size(); ++s) {
            move_item_to_dict("q" + std::to_string(s), ref_msm_data[s]);
            move_item_to_dict("psi" + std::to_string(s), ref_imt_data[s]);
        }
        return ref_return_dict;
    }

    [[ noreturn ]]
    void throw_value_error(char const *format, ...)
    {
        char buf[200];
        va_list args;
        va_start(args, format);
        std::vsnprintf(buf, sizeof(buf), format, args);
        va_end(args);
        throw std::domain_error(buf);
    }

    typedef PyObject *python_function(PyObject *, PyObject *, PyObject *);

    template <python_function FUNCTION>
    PyObject *parachute(PyObject *arg0, PyObject *args, PyObject *kwargs)
    {
        try {
            return FUNCTION(arg0, args, kwargs);
        }
        catch(std::domain_error const &e)
        {
            (void)PyErr_Format(PyExc_ValueError, "%s", e.what());
            return nullptr;
        }
        catch(std::exception const &e)
        {
            (void)PyErr_Format(PyExc_RuntimeError, "unexpected exception: %s", e.what());
            return nullptr;
        }
    }

#define WRAP_IN_PARACHUTE(function) ((PyCFunction)&parachute<function>)

    static int const MAX_S = MinkowskiAccumulator::MAX_S;

    static PyObject *wrap_imt_for_polygon(PyObject *, PyObject *args, PyObject *)
    {
        UniquePyPtr ref_vertices = nullptr;
        PyArrayObject *arr_vertices = nullptr;
        size_t N;

        {
            PyObject *arg1 = nullptr;
            if (!PyArg_ParseTuple(args, "O", &arg1, nullptr))
                return nullptr;
            ref_vertices = UniquePyPtr(
                PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));
            if (!ref_vertices)
                return nullptr;
        }

        arr_vertices = ref_vertices.reinterpret<PyArrayObject>();

        // make sure input has the correct dimensions
        if (PyArray_NDIM(arr_vertices) == 2) {
            npy_intp *dimensions = PyArray_DIMS(arr_vertices);
            N = dimensions[0];
            if (dimensions[1] != 2) {
                throw_value_error("data must be 2D, is %i", (int)dimensions[1]);
            } else if (N < 2) {
                throw_value_error("polygon must contain at least two vertices");
            }
        } else {
            throw_value_error("# dimensions must be 2, is %i", (int)PyArray_NDIM(arr_vertices));
        }

        std::vector<vec_t> vertices;
        vertices.reserve(N);

        for (size_t i = 0; i != N; ++i) {
            vertices.push_back({np_get_double(ref_vertices, i, 0),
                                np_get_double(ref_vertices, i, 1)});
        }

        auto return_data = MinkValReturnData(1, MAX_S);
        return_data.assign(0, papaya2::imt_polygon(vertices));
        return return_data.move_to_dict().release();
    }

#ifdef HAVE_CGAL
    static PyObject *wrap_imt_for_pointpattern(PyObject *, PyObject *args,
                                               PyObject *kwargs)
    {
        UniquePyPtr ref_seeds = nullptr;
        PyArrayObject *arr_seeds = nullptr;
        VoronoiDiagram vd;
        vd.periodic = false;
        int N;

        {
            PyObject *arg1 = nullptr;
            if (!PyArg_ParseTuple(args, "O", &arg1, nullptr))
                return nullptr;
            ref_seeds = UniquePyPtr(
                PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));
            if (!ref_seeds)
                return nullptr;
        }

        if (kwargs) for (auto key_and_value : Kwargs(kwargs))
        {
            // if box argument is given, assume we want a periodic Voronoi
            // tessellation
            if ("box" == key_and_value.first) {
                PyObject *boxarg = key_and_value.second;
                auto ref_box = UniquePyPtr(
                    PyArray_FROM_OTF(boxarg, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));
                PyArrayObject *arr_box = ref_box.reinterpret<PyArrayObject>();
                if (!arr_box) {
                    throw std::runtime_error("error converting box kwarg");
                }

                if (PyArray_NDIM(arr_box) == 1) {
                    npy_intp *dimensions = PyArray_DIMS(arr_box);
                    if (dimensions[0] == 2) {
                        vd.boxLx = np_get_double(ref_box, 0);
                        vd.boxLy = np_get_double(ref_box, 1);
                        vd.periodic = true;
                    } else if (dimensions[0] == 1) {
                        vd.boxLx = vd.boxLy = np_get_double(ref_box, 0);
                        vd.periodic = true;
                    } else {
                        throw_value_error("box kwarg must be length 1 oder length 2, is %i",
                            (int)dimensions[0]);
                    }
                } else if (PyArray_NDIM(arr_box) == 0) {
                    vd.boxLx = vd.boxLy = np_get_double(ref_box, 0);
                    vd.periodic = true;
                } else {
                    throw_value_error("box kwargs must be one-dimensional, is %i-dimensional",
                        (int)PyArray_NDIM(arr_box));
                }
            } else {
                throw_value_error("illegal keyword argument: %s", key_and_value.first.c_str());
            }
        }

        arr_seeds = ref_seeds.reinterpret<PyArrayObject>();
        // make sure input has the correct dimensions
        if (PyArray_NDIM(arr_seeds) == 2) {
            npy_intp *dimensions = PyArray_DIMS(arr_seeds);
            N = dimensions[0];
            if (dimensions[1] != 2) {
                throw_value_error("data must be 2D, is %i", (int)dimensions[1]);
            }
        } else {
            throw_value_error("# dimensions must be 2, is %i", (int)PyArray_NDIM(arr_seeds));
        }

        auto return_data = MinkValReturnData(N, MAX_S);

        // fill seed point coordinates into the Voronoi diagram
        for (int i = 0; i != N; ++i)
            vd.add_seed({np_get_double(ref_seeds, i, 0),
                         np_get_double(ref_seeds, i, 1)});

        // compute the Voronoi diagram
        vd.make_voronoi();

        {
            int label = 0;
            for (auto fit = vd.faces_begin(); fit != vd.faces_end(); ++fit) {
                VoronoiDiagram::point_t const seed = fit->dual()->point();

                if (seed != vd.seeds.at(label)) {
                    throw std::runtime_error("CGAL decided to reorder our seed points.");
                }

                // compute Minkowski valuations and copy them over
                if (!fit->is_unbounded()) {
                    return_data.assign(label, vd.minkval_for_cell(fit));
                }

                ++label;
                if (label == N)
                    break;
            }
        }

        return return_data.move_to_dict().release();
    }
#endif // HAVE_CGAL

    struct NumpyArrayPhoto
    {
        NumpyArrayPhoto(PyArrayObject *array)
            : array_(array)
        {
            npy_intp *dimensions = PyArray_DIMS(array);
            width_ = dimensions[0];
            height_ = dimensions[1];
        }

        unsigned long width_, height_;
        PyArrayObject *array_;

        double operator()(int i, int j) const
        {
            if (!(i >= 0 && i < width()))
                throw std::range_error("NumpyArrayPhoto i");
            if (!(j >= 0 && j < height()))
                throw std::range_error("NumpyArrayPhoto j");
            return *reinterpret_cast<double const *>(PyArray_GETPTR2(array_, i, j));
        }

        int width() const { return width_; }
        int height() const { return height_; }
        papaya2::vec_t origin() const { return {0., 0.}; }
        papaya2::vec_t upper_right() const
        {
            return {double(width()), double(height())};
        }
    };


    static PyObject *wrap_imt_for_image(PyObject *, PyObject *args,
                                        PyObject *kwargs)
    {
        UniquePyPtr ref_image = nullptr;
        double threshold = .5;
        double padding_value = 0;

        {
            PyObject *arg1 = nullptr;
            if (!PyArg_ParseTuple(args, "O", &arg1, nullptr))
                return nullptr;
            ref_image = UniquePyPtr(
                PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));
            if (!ref_image)
                return nullptr;
        }

        if (kwargs) for (auto key_and_value : Kwargs(kwargs))
        {
            if ("threshold" == key_and_value.first)
            {
                PyObject *threshold_arg = key_and_value.second;
                threshold = PyFloat_AsDouble(threshold_arg);

                if (threshold == -1. && PyErr_Occurred())
                    throw std::runtime_error ("error converting threshold kwarg");
            } else {
                throw_value_error("illegal keyword argument: %s", key_and_value.first.c_str());
            }
        }

        PyArrayObject *arr_image = ref_image.reinterpret<PyArrayObject>();
        if (PyArray_NDIM(arr_image) != 2) {
            throw_value_error("# dimensions of image must be 2, is %i",
                    (int)PyArray_NDIM(arr_image));
        }

        auto const original = NumpyArrayPhoto(arr_image);
        auto const padded = make_padded_view(original, padding_value);
        // FIXME: add periodic boundary conditions mode

        // get output arrays
        int const N = 1;

        auto return_data = MinkValReturnData(N, MAX_S);

        {
            MarchingSquaresFlags flags;
            auto minkval = imt_interpolated_marching_squares(padded, threshold, flags);
            int const label = 0;
            return_data.assign(label, minkval);
        }

        return return_data.move_to_dict().release();
    }

    static PyMethodDef mymethods[] = {
        {"imt_for_polygon", WRAP_IN_PARACHUTE(wrap_imt_for_polygon),
         METH_VARARGS,
         "imt_for_polygon(vertices)\n"
         "compute the Minkowski of the polygon bounded by the vertices. "
         "Vertices are assumed to be in counterclockwise order.\n"
         "Returns a dictionary mapping each metric to its value.\n"},
#ifdef HAVE_CGAL
        {"imt_for_pointpattern", WRAP_IN_PARACHUTE(wrap_imt_for_pointpattern),
         METH_VARARGS | METH_KEYWORDS,
         "imt_for_pointpattern(seeds)\n"
         "imt_for_pointpattern(seeds, box=[Lx, Ly])\n"
         "compute the Minkowski valuations of the Voronoi cells generated by "
         "the\n"
         "pointpattern.\n"
         "Returns a dictionary mapping each metric to a 1D NumPy array of "
         "length len(seeds).\n"
         "The elements of each array are the metrics corresponding to each "
         "seed point, in\n"
         "input order.  If any cell is unbounded, its metrics will be NaN.\n"},
#endif
        {"imt_for_image", WRAP_IN_PARACHUTE(wrap_imt_for_image),
         METH_VARARGS | METH_KEYWORDS,
         "imt_for_image(image)\n"
         "imt_for_image(seeds, threshold=value)\n"
         "Returns a dictionary that contains the metrics.\n"
         "The unit of scale is assumed to be 1 pixel.\n"},
        {nullptr, nullptr, 0, nullptr} // sentinel
    };
} // anonymous namespace

extern "C" {
#ifdef PYTHON_3
    static PyModuleDef module_def = {
        PyModuleDef_HEAD_INIT,
        "pypaya2",
        nullptr, // module docs
        -1,
        mymethods
    };

    PyMODINIT_FUNC PyInit_pypaya2()
    {
        import_array();
        return PyModule_Create(&module_def);
    }
#else // Python 2
    void initpypaya2()
    {
        (void)Py_InitModule("pypaya2", mymethods);
        import_array();
    }
#endif

} // extern "C"
