[![DOI](https://joss.theoj.org/papers/10.21105/joss.02538/status.svg)](https://doi.org/10.21105/joss.02538)
[![linux-test Status](https://github.com/morphometry/papaya2/workflows/linux-test/badge.svg)](https://github.com/morphometry/papaya2/actions)
[![mac-test Status](https://github.com/morphometry/papaya2/workflows/mac-test/badge.svg)](https://github.com/morphometry/papaya2/actions)

# Overview

papaya2 is a small header-only C++ library for computing irreducible 2D Minkowski tensors of image and polygonal data.
The library needs a C++ 11 compliant compiler, and no external dependencies.
More detailed information can be found at <https://morphometry.org/software/papaya2/>.

If you're using this work in published work, please cite

[Schaller et al., (2020). papaya2: 2D Irreducible Minkowski Tensor computation. Journal of Open Source Software, 5(54), 2538](https://doi.org/10.21105/joss.02538)

# Installation

papaya2 needs no installation, it is an header only library.


# Demos

papaya2 includes several demos which can be found in the demo folder,
see the [README](https://github.com/morphometry/papaya2/blob/master/demos/README.md) file.

An interactive demo can be found at <https://morphometry.org/morphometer>.
It is based on the JavaScript binding of papaya2.


# Tests

papaya2 inclues a test suite based on catch2. It can be run with

    cd test
    make


# Other language bindings

papaya2 includes bindings to Python, Matlab and JavaScript.
See READMEs in correspoding folders for more details or <https://morphometry.org/software/papaya2/>.


# Contributing

If you have some contribution to the software, please write an email to <info@morphometry.org>.
Bugs can be filed on [github](https://github.com/morphometry/papaya2/issues).
