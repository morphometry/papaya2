# Usage

# Installation

Install the header files for Python 3, numpy and CGAL.
(On Debian-based systems: `apt-get install python3-dev python3-numpy cgal-dev`.)

To compile type `make`.

You can put the resulting `pypaya2.so` in the same directory as your analysis
Python scripts, or on your PYTHONPATH.

# Usage

Computing the Minkowski Tensors of a polygon:

    import pypaya2
    vertices = [[2., 3.], [1., 1.], [-3., 4.], [-2., 0.], [2., 0.]
    minkval = pypaya2.imt_for_polygon(vertices)

Performing Voronoi-Minkowski analysis of a point pattern:

    import pypaya2
    seeds = [[2., 3.], [1., 1.], [-3., 4.], [-2., 0.], [2., 0.]
    minkval = pypaya2.imt_for_pointpattern(seeds)

There is currently no support for image analysis.  Please look to the
demos folder at the toplevel if you need that.

## Compiling for Python 2

If you still have Python 2, go to features.mk and add the line `PYTHON_VERSION = 2`.
You will need `python2-dev` and `python-numpy` packages.

## Running the tests

Run `python3 test.py`.

## Common errors

If you get the error message

    Could not find the numpy C headers. Is the numpy Python module installed properly?

You do not have `numpy` installed.  On Debian-based systems: `sudo apt-get install python3-dev python3-numpy`.
On Macintosh with Homebrew: `brew install numpy`.  On Macintosh with Conda: `conda install -c anaconda numpy`.
