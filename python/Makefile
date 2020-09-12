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
    PYTHON_CXXFLAGS += `python3 ../config/detect-numpy.py`
else
    PYTHON_CXXFLAGS += `python-config --includes`
    PYTHON_CXXFLAGS += `python ../config/detect-numpy.py`
endif

# override defaults here
-include features.mk

CXXFLAGS += -I../include $(PYTHON_CXXFLAGS)
LDFLAGS += $(PYTHON_LDFLAGS) $(CGAL_LDFLAGS)

PRODUCTS = \
    pypaya2.$(SO_EXTENSION) \

default: pypaya2

include ../config/shared.mk
include ../config/detect-cgal.mk

all: pypaya2

pypaya2: pypaya2.$(SO_EXTENSION)
	@true

pypaya2.$(SO_EXTENSION): pypaya2.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $< $(LDFLAGS)
