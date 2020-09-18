// point pattern analysis: compute Voronoi diagram for a set
// of seed points, and the IMTs of the Voronoi cells.
// optionally, with periodic boundary conditions
// 2019 Jenny Wagner <j.wagner@uni-heidelberg.de>
// 2019 Sebastian Kapfer <sebastian.kapfer@fau.de>
// 2019 Fabian Schaller <physik@fabian-schaller.de>
// FIXME this does not have tests yet
//
//
// Compiling
//
// Compiling this in MATLAB is a little tricky because we use the CGAL
// library.
//
// If you have CGAL with a version < 5.0
//
// mex -v -I"../include" -lCGAL -lCGAL_Core -lgmp imt_for_pointpattern.cpp
//
// if you have CGAL with a version >= 5.0
//
// mex -v CXXFLAGS='-std=c++14' -I"../include" -lgmp imt_for_pointpattern.cpp
//
// If in doubt, try both.
//
// If the CGAL library is not installed in one of the paths searched
// by your C++ compiler, you may need to add paths to CGAL, like this:
// ... -I<path_to_CGAL> -L<path_to_CGAL_libs> ...
//
//
// Usage
//
// function call for periodic boundaries:
// output_array = imt_for_pointpattern(input_array, [box_side_x, box_side_y])
//
// function call for non-periodic boundaries:
// output_array = imt_for_pointpattern(input_array)
//
// input_array is a  n x 2 matrix containing the seed point coordinates (x,y)
// output_array is a n x 9 matrix of double values with columns
//
// col1: seedx
// col2: seedy
// col3: area
// col4: perimeter
// col5: q2
// col6: q3
// col7: q4
// col8: q5
// col9: q6

#include <papaya2/voronoi.hpp>
#include <gmpxx.h>
#include <mex.hpp>
#include <mexAdapter.hpp>

using namespace papaya2;
using papaya2::VoronoiDiagram;
using point_t = VoronoiDiagram::point_t;

// the interface function to MATLAB
struct MexFunction : matlab::mex::Function
{
    using ArgumentList = matlab::mex::ArgumentList;

    void operator()(ArgumentList outputs, ArgumentList inputs)
    {
        VoronoiDiagram vd;

        matlab::data::Array const seeds = std::move(inputs[0]);
        auto const dimensions = seeds.getDimensions();
        std::cerr << "read in dimensions of input array: " << dimensions[0]
                  << " " << dimensions[1] << std::endl;

        for (int i = 0; i < dimensions[0]; ++i)
            vd.add_seed({seeds[i][0], seeds[i][1]});

        if (inputs.size() > 1) {
            matlab::data::Array const box = std::move(inputs[1]);
            vd.boxLx = box[0];
            vd.boxLy = box[1];
            vd.periodic = true;
        }

        matlab::data::ArrayFactory f;
        matlab::data::TypedArray<double> outarray =
            f.createArray<double>({dimensions[0], 9});

        if (vd.periodic) {
            point_t offender;
            if (!vd.seeds_in_box(&offender)) {
                std::cerr << "Warning: There are points outside of the box: "
                          << offender.x() << ", " << offender.y() << "\n";
            }
        }

        vd.make_voronoi();

        int label = 0;
        for (auto fit = vd.faces_begin(); fit != vd.faces_end(); ++fit) {
            point_t const seed = fit->dual()->point();

            if (seed != vd.seeds.at(label)) {
                throw std::runtime_error("CGAL decided to reorder our seed points.");
            }

            if (!fit->is_unbounded()) {
                auto minkval = vd.minkval_for_cell(fit);

                outarray[label][0] = seed.x();
                outarray[label][1] = seed.y();
                outarray[label][2] = minkval.area();
                outarray[label][3] = minkval.perimeter();

                for (auto s : {2, 3, 4, 5, 6}) {
                    outarray[label][s + 2] = minkval.msm(s);
                }
            } else {
                outarray[label][0] = seed.x();
                outarray[label][1] = seed.y();
                for (auto s : {2, 3, 4, 5, 6}) {
                    outarray[label][s + 2] = 0.0 / 0.0;
                }
            }

            ++label;
            if (label == (int)vd.seeds.size())
                break;
        }

        outputs[0] = outarray;
    }
};
