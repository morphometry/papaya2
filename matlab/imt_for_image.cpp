// a MATLAB extension analyzing a 2D image using irreducible
// Minkowski tensors using the preliminary Papaya2 library.
// for an example, see imganalysis.m
// 2019 Sebastian Kapfer <sebastian.kapfer@fau.de>
// 2019 Jenny Wagner <j.wagner@uni-heidelberg.de>
//
// compile in matlab with the command
// mex -v -I ../include imt_for_image.cpp
//
// function call:
// output_array = imt_for_image(input_image, thresholds)
//
// output_array is a numt x 17 matrix of double values with columns
// col1: threshold
// col2: area
// col3: perimeter
// col4: q_2
// col5: arg(psi_2)
// col6: q_3
// col7: arg(psi_3) ...
// col17: arg(psi_8)

#include <papaya2.hpp>
#include <iostream>
#include <mex.hpp>
#include <mexAdapter.hpp>

using namespace papaya2;

// structure to implement data transfer from Matlab to C++.
// input argument is a Matlab data structure containing an array of doubles.
struct MexPhoto
{
    matlab::data::Array data_;
    int width_, height_;

    MexPhoto(matlab::data::Array &&data) : data_(data)
    {
        // find dimensions of input
        auto const dimensions = data_.getDimensions();
        if (dimensions.size() != 2) {
            throw std::runtime_error("input image must be 2D");
        }
        width_ = dimensions[1];
        height_ = dimensions[0];
    }

    double operator()(int i, int j) const
    {
        if (i < 0 || i >= width_ || j < 0 || j >= height_)
            throw std::range_error("invalid pixel index in MexPhoto");
        return data_[height_ - 1 - j][i];
    }

    int width() const { return width_; }
    int height() const { return height_; }
    vec_t origin() const { return { 0., 0. }; }
    vec_t upper_right() const { return { double(width_), double(height_) }; }
};

// the interface function to MATLAB
struct MexFunction : matlab::mex::Function
{
    // define function for banana-MATLAB-call
    void operator()(matlab::mex::ArgumentList outputs,
                    matlab::mex::ArgumentList inputs)
    {

        // process arguments of call
        MexPhoto image(std::move(inputs[0]));

        std::vector<double> thresholds;

        {
            matlab::data::Array arg2 = std::move(inputs[1]);

            auto const dimensions = arg2.getDimensions();
            if (dimensions.size() != 2 || dimensions[0] != 1) {
                throw std::runtime_error("thresholds must be 1D array");
            }

            for (int i = 0; i != dimensions[1]; ++i)
                thresholds.push_back(arg2[i]);
        }

        // assemble output array for matlab
        int n = 0;

        matlab::data::ArrayFactory f;

        matlab::data::TypedArray<double> outarray =
            f.createArray<double>({thresholds.size(), 17});

        for (double thresh : thresholds) {
            std::cerr << "IMT analysis for threshold[" << n << "] = " << thresh
                      << "\n";
            auto imt =
                papaya2::imt_interpolated_marching_squares(image, thresh);

            outarray[n][0] = thresh;
            outarray[n][1] = imt.area();
            outarray[n][2] = imt.perimeter();
            outarray[n][3] = imt.msm(2);
            outarray[n][4] = std::arg(imt.imt(2));
            outarray[n][5] = imt.msm(3);
            outarray[n][6] = std::arg(imt.imt(3));
            outarray[n][7] = imt.msm(4);
            outarray[n][8] = std::arg(imt.imt(4));
            outarray[n][9] = imt.msm(5);
            outarray[n][10] = std::arg(imt.imt(5));
            outarray[n][11] = imt.msm(6);
            outarray[n][12] = std::arg(imt.imt(6));
            outarray[n][13] = imt.msm(7);
            outarray[n][14] = std::arg(imt.imt(7));
            outarray[n][15] = imt.msm(8);
            outarray[n][16] = std::arg(imt.imt(8));

            ++n;
        }

        outputs[0] = outarray;
    }
};
