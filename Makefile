# defaults
CXXFLAGS += -g -O3
CXXFLAGS += -Wall -Werror
# enable FITS support being able to load FITS containers (in banana)
FITS_SUPPORT = 1
FITS_LDFLAGS = -lCCfits -lcfitsio
# enable CGAL support for Voronoi diagrams (in ppanalysis and in pypaya2)
CGAL_SUPPORT = 1
CGAL_LDFLAGS = -lCGAL -LCGAL_Core -lgmp
# set to "dll" for Windows
SO_EXTENSION = so
CLANGFORMAT = clang-format

# modify defaults here if necessary
-include features.mk

# FITS support
ifeq ($(FITS_SUPPORT), 1)
    CXXFLAGS += -DHAVE_CCFITS
    LDFLAGS += $(FITS_LDFLAGS)
endif

# have CGAL for Voronoi diagrams?
ifeq ($(CGAL_SUPPORT), 1)
    CXXFLAGS += -DHAVE_CGAL
    LDFLAGS += $(CGAL_LDFLAGS)
endif

CXXFLAGS += -std=c++11 -fPIC

MAKEFILES = \
    Makefile \
    features.mk \

TEST_OBJECTS = \
    test_image.o \
    test_tools.o \
    test_imt.o \

BINARIES = \
    banana \
    runtests \
    sersic \
    imganalysis \
    ppanalysis \
    pypaya2.$(SO_EXTENSION) \

default: all runtests
	./runtests

all: .ts.mk.hpp imganalysis ppanalysis pypaya2.$(SO_EXTENSION)

# hack to make any code recompile if Makefile changes
.ts.mk.hpp: $(MAKEFILES)
	@touch $@

# automagic dependencies
depend: *.cpp *.hpp .ts.mk.hpp
	$(CXX) $(CXXFLAGS) -MM *.cpp >.depend
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

pypaya2: pypaya2.$(SO_EXTENSION)
	@true

pypaya2.$(SO_EXTENSION): pypaya2.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $< $(LDFLAGS)

# special rule, needs Python headers
pypaya2.o: pypaya2.cpp
	$(CXX) $(CXXFLAGS) `python-config --includes` -c -o $@ $<

test: runtests
	./runtests

pretty:
	$(CLANGFORMAT) -i -style=file `ls *pp | grep -vE '(catch.hpp|picopng.hpp)'`

clean:
	rm -f $(BINARIES) *.o .ts.mk.hpp

.PHONY: default all test clean pretty depend
