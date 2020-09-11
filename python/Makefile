# defaults
CXXFLAGS += -g -O3
CXXFLAGS += -Wall -Werror
CXXFLAGS += -fPIC
# set to "dll" for Windows
SO_EXTENSION = so
# compile pypaya2 for Python 3
PYTHON_VERSION = 3

ifeq ($(PYTHON_VERSION), 3)
    PYTHON_CXXFLAGS += -DPYTHON_3
    PYTHON_CXXFLAGS += `python3-config --includes`
    PYTHON_CXXFLAGS += -I/usr/lib/python3/dist-packages/numpy/core/include
else
    PYTHON_CXXFLAGS += `python-config --includes`
endif

include ../cgal.mk

# override defaults here
-include features.mk

CXXFLAGS += -I../include $(PYTHON_CXXFLAGS) $(CGAL_CXXFLAGS)
LDFLAGS += $(PYTHON_LDFLAGS) $(CGAL_LDFLAGS)

PRODUCTS = \
    pypaya2.$(SO_EXTENSION) \

default: .ts.mk.hpp pypaya2.$(SO_EXTENSION)

include ../config/shared.mk

pypaya2: pypaya2.$(SO_EXTENSION)
	@true

pypaya2.$(SO_EXTENSION): pypaya2.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $< $(LDFLAGS)
