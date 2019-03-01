# defaults
CXXFLAGS += -g -O3
CXXFLAGS += -Wall -Werror
FITS_SUPPORT = 0
CGAL_SUPPORT = 1
BUILD_PYTHON_MODULE = 1

# modify defaults here if necessary
-include features.mk

# FITS support
ifeq ($(FITS_SUPPORT), 1)
    CXXFLAGS += -DHAVE_CCFITS
    LDFLAGS += -lCCfits -lcfitsio
endif

# have CGAL for Voronoi diagrams?
ifeq ($(CGAL_SUPPORT), 1)
    CXXFLAGS += -DHAVE_CGAL
    LDFLAGS += -lCGAL -LCGAL_Core -lgmp
endif

# do we want to build the Python2 module?
ifeq ($(BUILD_PYTHON_MODULE), 1)
    CXXFLAGS += -I/usr/include/python2.7
endif

CXXFLAGS += -std=c++11 -fPIC

MAKEFILES = \
    Makefile \
    features.mk \

TEST_SOURCES = \
    test_image.cpp \
    test_tools.cpp \
    test_imt.cpp \

TEST_OBJECTS = $(TEST_SOURCES:%.cpp=%.o)

ALL_SOURCES = \
    $(TEST_SOURCES) \
    runtests.cpp \
    sersic.cpp \
    ppanalysis.cpp \
    pypaya2.cpp \

BINARIES = \
    banana \
    runtests \
    sersic \
    ppanalysis \
    pypaya2.so \

default: all runtests
	./runtests

all: .ts.mk.hpp .depend imganalysis ppanalysis pypaya2.so

# hack to make any code recompile if Makefile changes
.ts.mk.hpp: $(MAKEFILES)
	@touch $@

# automagic dependencies
.depend: $(ALL_SOURCES) *.hpp .ts.mk.hpp
	$(CXX) $(CXXFLAGS) -MM $(ALL_SOURCES) >.depend
-include .depend

# create features.mk if absent
features.mk:
	@touch $@

imganalysis: imganalysis.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

banana: banana.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

ppanalysis: ppanalysis.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

runtests: $(TEST_OBJECTS) runtests.o
	$(CXX) $(CXXFLAGS) -o $@ $(TEST_OBJECTS) runtests.o $(LDFLAGS)

sersic: sersic.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

pypaya2.so: pypaya2.o
	$(CXX) -shared -o $@ $< $(LDFLAGS)

test: runtests
	./runtests

pretty:
	$(CLANGFORMAT) -i -style=file `ls *pp | grep -vE '(catch.hpp|picopng.hpp)'`

clean:
	rm -f $(BINARIES) *.o .ts.mk.hpp

.PHONY: default all test clean
