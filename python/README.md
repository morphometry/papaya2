# Usage

# Installation

1. Install the header files for Python 3, numpy and CGAL (see the end of this document).
On Debian-based systems: `sudo apt-get install libcgal-dev python3-dev python3-numpy`.

2. To compile type `make`.

   If you get the error message

        Could not find the numpy C headers. Is the numpy Python module installed properly?

   you do not have the [`numpy`](https://numpy.org/) Python package installed.

3. You can put the resulting `pypaya2.so` in the same directory as your analysis
Python scripts, or on your PYTHONPATH.

# Usage

Computing the Minkowski Tensors of a polygon:

    import pypaya2
    vertices = [[2., 3.], [1., 1.], [-3., 4.], [-2., 0.], [2., 0.]
    minkval = pypaya2.imt_for_polygon(vertices)
    print(minkval['psi2'])

Performing Voronoi-Minkowski analysis of a point pattern:

    import pypaya2
    seeds = [[2., 3.], [1., 1.], [-3., 4.], [-2., 0.], [2., 0.]
    minkval = pypaya2.imt_for_pointpattern(seeds)
    print(minkval['psi2'])

There is currently no support for image analysis.  Please look to the
demos folder at the toplevel if you need that.

## Running the tests

Run `make test`.

## Installing dependencies

On Debian/Ubuntu-based systems:

    sudo apt-get install libcgal-dev python3-dev python3-numpy

On MacOS with [Homebrew](https://docs.brew.sh/)

    brew install boost cgal numpy

**Once you have installed any missing dependencies, type `make clean` to have the Makefile detect them.**

## Compiling for Python 2

If you still have Python 2, edit the file `features.mk` and add the line `PYTHON_VERSION = 2`.
You will need `python2-dev` and `python-numpy` packages.
