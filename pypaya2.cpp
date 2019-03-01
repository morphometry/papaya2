// Python/Numpy interface.
// we're not using Boost::Python here as we export very few
// functions and Boost::Python::Numpy is not in Debian.
// 2019 Sebastian Kapfer <sebastian.kapfer@fau.de>
#include "papaya2.hpp"
#ifdef HAVE_CGAL
// without CGAL, there is no Voronoi.
#include "voronoi.hpp"
#endif
#include "python.hpp"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

using namespace papaya2;
using namespace papaya2::python;

#ifdef HAVE_CGAL
static void np_set_double(PyArrayObject *arr, int i, double value)
{
    *reinterpret_cast<double *>(PyArray_GETPTR1(arr, i)) = value;
}

static void np_set_complex(PyArrayObject *arr, int i, complex_t value)
{
    *reinterpret_cast<complex_t *>(PyArray_GETPTR1(arr, i)) = value;
}

static double np_get_double(PyArrayObject *arr, int i)
{
    return *reinterpret_cast<double const *>(PyArray_GETPTR1(arr, i));
}

static double np_get_double(PyArrayObject *arr, int i, int j)
{
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

enum
{
    METR_PERI,
    METR_Q2,
    METR_Q3,
    METR_Q4,
    METR_Q5,
    METR_Q6,
    METR_Q7,
    METR_Q8,
    METR_PSI2,
    METR_PSI6,
    NUM_METR
};

#endif // HAVE_CGAL

extern "C"
{

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

        if (kwargs) {
            // FIXME flag illegal kwargs

            // if box argument is given, assume we want a periodic Voronoi
            // tessellation
            if (PyObject *boxarg = PyDict_GetItemString(kwargs, "box")) {
                auto ref_box = UniquePyPtr(
                    PyArray_FROM_OTF(boxarg, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));
                PyArrayObject *arr_box = ref_box.reinterpret<PyArrayObject>();
                if (!arr_box) {
                    (void)PyErr_Format(PyExc_RuntimeError,
                                       "error converting box kwarg");
                    return nullptr;
                }

                if (PyArray_NDIM(arr_box) == 1) {
                    npy_intp *dimensions = PyArray_DIMS(arr_box);
                    if (dimensions[0] == 2) {
                        vd.boxLx = np_get_double(arr_box, 0);
                        vd.boxLy = np_get_double(arr_box, 1);
                        vd.periodic = true;
                    } else if (dimensions[0] == 1) {
                        vd.boxLx = vd.boxLy = np_get_double(arr_box, 0);
                        vd.periodic = true;
                    } else {
                        (void)PyErr_Format(
                            PyExc_ValueError,
                            "box kwarg must be length 1 oder length 2, is %i",
                            (int)dimensions[0]);
                        return nullptr;
                    }
                } else if (PyArray_NDIM(arr_box) == 0) {
                    vd.boxLx = vd.boxLy = np_get_double(arr_box, 0);
                    vd.periodic = true;
                } else {
                    (void)PyErr_Format(
                        PyExc_ValueError,
                        "box kwargs must be one-dimensional, is %i-dimensional",
                        (int)PyArray_NDIM(arr_box));
                    return nullptr;
                }
            }
        }

        arr_seeds = ref_seeds.reinterpret<PyArrayObject>();
        // make sure input has the correct dimensions
        if (PyArray_NDIM(arr_seeds) == 2) {
            npy_intp *dimensions = PyArray_DIMS(arr_seeds);
            N = dimensions[0];
            if (dimensions[1] == 2) {
                // OK
            } else {
                (void)PyErr_Format(PyExc_ValueError, "data must be 2D, is %i",
                                   (int)dimensions[1]);
                return nullptr;
            }
        } else {
            // PyErr_Format always returns 0 by documentation
            (void)PyErr_Format(PyExc_ValueError,
                               "# dimensions must be 2, is %i",
                               (int)PyArray_NDIM(arr_seeds));
            return nullptr;
        }

        // get output arrays
        auto ref_return_data = std::vector<UniquePyPtr>(NUM_METR);
        for (auto m : {METR_PERI, METR_Q2, METR_Q3, METR_Q4, METR_Q5, METR_Q6,
                       METR_Q7, METR_Q8})
            ref_return_data[m] = np_make_new_vector(N, NPY_DOUBLE);
        for (auto m : {METR_PSI2, METR_PSI6})
            ref_return_data[m] = np_make_new_vector(N, NPY_CDOUBLE);

        // fill seed point coordinates into the Voronoi diagram
        for (int i = 0; i != N; ++i)
            vd.add_seed({np_get_double(arr_seeds, i, 0),
                         np_get_double(arr_seeds, i, 1)});

        // compute the Voronoi diagram
        vd.make_voronoi();

        {
            int label = 0;
            for (auto fit = vd.faces_begin(); fit != vd.faces_end(); ++fit) {
                VoronoiDiagram::point_t const seed = fit->dual()->point();

                if (seed != vd.seeds.at(label)) {
                    (void)PyErr_Format(
                        PyExc_RuntimeError,
                        "CGAL decided to reorder our seed points.");
                    return nullptr;
                }

                // compute Minkowski valuations and copy them over
                if (!fit->is_unbounded()) {
                    auto const minkval = vd.minkval_for_cell(fit);
                    np_set_double(
                        ref_return_data[METR_PERI].reinterpret<PyArrayObject>(),
                        label, minkval.perimeter());
                    np_set_double(
                        ref_return_data[METR_Q2].reinterpret<PyArrayObject>(),
                        label, minkval.msm(2));
                    np_set_double(
                        ref_return_data[METR_Q3].reinterpret<PyArrayObject>(),
                        label, minkval.msm(3));
                    np_set_double(
                        ref_return_data[METR_Q4].reinterpret<PyArrayObject>(),
                        label, minkval.msm(4));
                    np_set_double(
                        ref_return_data[METR_Q5].reinterpret<PyArrayObject>(),
                        label, minkval.msm(5));
                    np_set_double(
                        ref_return_data[METR_Q6].reinterpret<PyArrayObject>(),
                        label, minkval.msm(6));
                    np_set_double(
                        ref_return_data[METR_Q7].reinterpret<PyArrayObject>(),
                        label, minkval.msm(7));
                    np_set_double(
                        ref_return_data[METR_Q8].reinterpret<PyArrayObject>(),
                        label, minkval.msm(8));
                    np_set_complex(
                        ref_return_data[METR_PSI2].reinterpret<PyArrayObject>(),
                        label, minkval.imt(2));
                    np_set_complex(
                        ref_return_data[METR_PSI6].reinterpret<PyArrayObject>(),
                        label, minkval.imt(6));
                }

                ++label;
                if (label == N)
                    break;
            }
        }

        // allocate the return dict and populate it
        auto ref_return_dict = UniquePyPtr(PyDict_New());
        std::array<char const *, NUM_METR> const metr_name = {
            "P", "q2", "q3", "q4", "q5", "q6", "q7", "q8", "psi2", "psi6"};
        for (int j = 0; j != NUM_METR; ++j)
            PyDict_SetItemString(ref_return_dict.get(), metr_name[j],
                                 ref_return_data[j].release());

        return ref_return_dict.release();
    }
#endif // HAVE_CGAL

    static PyMethodDef mymethods[] = {
#ifdef HAVE_CGAL
        {"imt_for_pointpattern", (PyCFunction)wrap_imt_for_pointpattern,
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
        {nullptr, nullptr, 0, nullptr} // sentinel
    };

    void initpypaya2()
    {
        (void)Py_InitModule("pypaya2", mymethods);
        import_array();
    }

} // extern "C"
