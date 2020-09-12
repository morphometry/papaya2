# defaults
CXXFLAGS += -g -O3 -std=c++11
CXXFLAGS += -Wall -Werror
CXXFLAGS += -fPIC
# set to "dll" for Windows
SO_EXTENSION = so
# compile pypaya2 for Python 3
PYTHON_VERSION = 3

# override defaults here
-include features.mk

ifeq ($(PYTHON_VERSION), 3)
    PYTHON = python3
    PYTHON_CXXFLAGS += -DPYTHON_3
    PYTHON_CXXFLAGS += `python3-config --includes`
    PYTHON_CXXFLAGS += `python3 ../config/detect-numpy.py`
    PYTHON_LDFLAGS += `python3-config --ldflags`
else
    PYTHON = python2
    PYTHON_CXXFLAGS += `python-config --includes`
    PYTHON_CXXFLAGS += `python ../config/detect-numpy.py`
    PYTHON_LDFLAGS += `python-config --ldflags`
endif

CXXFLAGS += -I../include $(PYTHON_CXXFLAGS)
LDFLAGS += $(PYTHON_LDFLAGS) $(CGAL_LDFLAGS)

DYNLIB_FLAGS = -shared
ifeq ("$(shell uname)", "Darwin")
    # unresolved symbols are Python's and will get resolved when loading
    DYNLIB_FLAGS += -undefined dynamic_lookup
endif

PRODUCTS = \
    pypaya2.$(SO_EXTENSION) \

default: pypaya2

include ../config/shared.mk
include ../config/detect-cgal.mk

all: pypaya2

pypaya2: pypaya2.$(SO_EXTENSION)
	@true

pypaya2.$(SO_EXTENSION): pypaya2.o
	$(CXX) $(CXXFLAGS) $(DYNLIB_FLAGS) -o $@ $< $(LDFLAGS)

non_cgal_test: pypaya2
	cp -f pypaya2.$(SO_EXTENSION) test
	(cd test; $(PYTHON) poly_test.py)

test: non_cgal_test
	(cd test; $(PYTHON) pp_test.py)
