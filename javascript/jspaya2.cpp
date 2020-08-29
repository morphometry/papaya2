  // 2020 Fabian Schaller <physik@fabian-schaller.de>

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include "papaya2.hpp"

using namespace emscripten;
using namespace papaya2;


class Picture {
public:
  Picture(){}

  using data_t = double;

  void initialize(int a,int b,double value){
    sx = a;
    sy = b;
    grayscaleImg.resize(sx*sy);
    for (unsigned int i = 0; i < sx*sy; ++i){
      grayscaleImg.at(i) = value;
    }
  }
  double get(int x, int y) const {
    unsigned int pos = (x + y*sx);
    return grayscaleImg.at(pos);
  }
  double operator()(int x, int y) const
  {
    return get(x, y);
  }
  int width() const
  {
    return sx;
  }
  int height() const
  {
    return sy;
  }
  vec_t origin() const
  {
    return ORIGIN;
  }
  vec_t upper_right() const
  {
    return {static_cast<double>(sx),static_cast<double>(sy)};
  }
  double pixel_width() const { return 1.; }
  double pixel_height() const { return 1.; }
  void readFromCanvasContextImgData(val array){
    auto v = vecFromJSArray<double>(array);
    for (int i = 0; i < v.size(); i=i+4){
      grayscaleImg.at(i/4)=v.at(i);
    }
  }

  MinkowskiAccumulator ni_minkowskis (double seg_threshold, bool analyze_white = true) const
  {
    if (analyze_white == true)
      return imt_regular_marching_squares(*this, seg_threshold);
    else
      return imt_regular_marching_squares(*this, seg_threshold, ANALYZE_BLACK | CONNECT_BLACK);
  }

  MinkowskiAccumulator minkowskis (double seg_threshold, bool analyze_white = true) const
  {
    if (analyze_white == true)
      return imt_interpolated_marching_squares(*this, seg_threshold);
    else
      return imt_interpolated_marching_squares(*this, seg_threshold, ANALYZE_BLACK | CONNECT_BLACK);
  }

private:
  std::vector<double> grayscaleImg;
  int sx;
  int sy;
};



MinkowskiAccumulator imt_single_poly(val vertices){
  auto v_points = vecFromJSArray<vec_t>(vertices);
  MinkowskiAccumulator macc;
  add_polygon_area(&macc, v_points);
  add_polygon_contour(&macc, v_points);
  return macc;
}


// Binding code
EMSCRIPTEN_BINDINGS(my_class_example) {
  value_array<vec_t>("vec_t")
    .element(emscripten::index<0>())
    .element(emscripten::index<1>())
    ;
  class_<MinkowskiAccumulator>("MinkowskiAccumulator")
    .constructor()
    .function("area", &MinkowskiAccumulator::area)
    .function("peri", &MinkowskiAccumulator::perimeter)
    .function("msm", &MinkowskiAccumulator::msm)
    .function("imt", &MinkowskiAccumulator::imt)
    .function("beta102", &MinkowskiAccumulator::beta102)
    .function("isoper" , &MinkowskiAccumulator::isoper)
    ;
  class_<complex_t>("Complex")
    .constructor()
    .function("real", static_cast< double (complex_t::*)() const >(&complex_t::real))
    .function("imag", static_cast< double (complex_t::*)() const >(&complex_t::imag))
    ;
  class_<Picture>("Picture")
    .constructor()
    .function("initialize", &Picture::initialize)
    .function("get", &Picture::get)
    .function("readFromCanvasContextImgData", &Picture::readFromCanvasContextImgData)
    .function("noninterpolated_minkowskis", &Picture::ni_minkowskis)
    .function("interpolated_minkowskis", &Picture::minkowskis)
    ;
  function("imt_single_poly", &imt_single_poly);
  register_vector<int>("vector<int>");
  //vec_t from papaya2 tools.hpp
  register_vector<vec_t>("vector<vec_t>");
  register_map<int,MinkowskiAccumulator>("MapIntMinkAcc");
  register_map<int,std::string>("map<int, string>");
}


