#pragma once

// example of how to write a class suitable for the PHOTO parameter
// of the template functions in namespace papaya2.
// this example loads the red color channel from a PNG file.
// it is used in imganalysis.cpp.

#include "picopng.hpp"

namespace {

using bytevec_t = std::vector<unsigned char>;
using string = std::string;

static void load_file(bytevec_t *buffer, const string &filename)
{
    std::ifstream file(filename.c_str(),
                       std::ios::in | std::ios::binary | std::ios::ate);

    // get filesize
    std::streamsize size = 0;
    if (file.seekg(0, std::ios::end).good())
        size = file.tellg();
    if (file.seekg(0, std::ios::beg).good())
        size -= file.tellg();

    // read contents of the file into the vector
    if (size > 0) {
        buffer->resize((size_t)size);
        file.read((char *)buffer->data(), size);
    } else
        buffer->clear();
}

struct PingFile
{
    unsigned long width_, height_;
    bytevec_t buffer_;

    enum
    {
        RED,
        GREEN,
        BLUE,
        ALPHA
    };

    PingFile(const string &infilename)
    {
        bytevec_t png_file_data;
        load_file(&png_file_data, infilename);
        if (!png_file_data.size())
            throw std::runtime_error("file is empty");
        int error = decodePNG(buffer_, width_, height_, png_file_data.data(),
                              png_file_data.size(), true);
        if (error != 0)
            throw std::runtime_error("picopng error " + std::to_string(error));
        if (buffer_.size() != width_ * height_ * 4)
            throw std::runtime_error("picopng didn't decode right");
    }

    int operator()(int i, int j) const
    {
        return const_cast<PingFile *>(this)->at(i, j, RED);
    }

    unsigned char &at(int i, int j, int c)
    {
        if (!(i >= 0 && i < width()))
            throw std::range_error("PingFile i");
        if (!(j >= 0 && j < height()))
            throw std::range_error("PingFile j");
        if (!(c >= 0 && c < 4))
            throw std::range_error("PingFile c");
        j = height() - 1 - j;
        int index = (i + j * width()) * 4 + c;
        return buffer_.at(index);
    }

    int width() const { return width_; }
    int height() const { return height_; }
    papaya2::vec_t origin() const { return {0., 0.}; }
    papaya2::vec_t upper_right() const
    {
        return {double(width()), double(height())};
    }
};

} // namespace
