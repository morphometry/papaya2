# Demos

There are four demos explaining how to use the `papaya2` library in this directory.
They are meant to be modified and adapted to your needs.
It should be pretty easy to get them to compile on Linux/BSD/Mac systems.

More detailed information can be found at <https://morphometry.org/software/papaya2/>.

## `imganalysis`

`imganalysis` demonstrates how to analyze pixel data.  Type `make imganalysis` to build it;
it needs no external dependencies.

`imganalysis` reads images in PNG format.

    ./imganalysis in example_inputs/GRF_matern_C5.png out data.txt

See `./imganalysis help` and <https://morphometry.org/software/papaya2/> for further information.

## `ppanalysis`

`ppanalysis` analyzes point patterns.  For building it the CGAL library is required,
see the end of this document.  Type `make ppanalysis` to compile.

An example file can be run by executing

    ./ppanalysis in example_inputs/granular-cryst-cluster.txt out out.txt boxL 500    

See `./ppanalysis help` and <https://morphometry.org/software/papaya2/> for further information.

## `banana`

`banana` analyzes astrophysics data in FITS format.  The `CCfits` library is required to build,
see below.

    ./banana in example_inputs/test_SIE_detail2.fits out test_SIE_detail2_out.dat mint 0.0003 maxt 0.003 numt 100

See `./banana help` and <https://morphometry.org/software/papaya2/> for further information.

## `sersic`

`sersic` is another example analyzing pixel data.  It samples 
[Sérsic density profiles](https://en.wikipedia.org/wiki/Sersic_profile)
and computes Minkowski Tensors of their excursion sets.

    ./sersic scan_angle threshold 1.5 aspect 0.3 resolution 300 interpolated_marching_squares

The `sersic` example does not require external dependencies.

## Installing dependencies

On Debian/Ubuntu-based systems:

    sudo apt-get install libcgal-dev libccfits-dev

On MacOS with [Homebrew](https://docs.brew.sh/)

    brew install boost cgal ccfits

Once you have installed any missing dependencies, type `make clean` to have the Makefile detect them.
