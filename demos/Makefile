# defaults
CXXFLAGS += -g -O3
CXXFLAGS += -Wall -Werror
FITS_LDFLAGS = -lCCfits -lcfitsio
CXXFLAGS += -I../include

include ../cgal.mk

# modify defaults here if necessary
-include features.mk

PRODUCTS = \
    banana \
    sersic \
    imganalysis \
    ppanalysis \

default: .ts.mk.hpp imganalysis ppanalysis

include ../config/shared.mk

all: default banana sersic

imganalysis: imganalysis.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

# special rule, to be able to load FITS containers
banana: banana.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(FITS_LDFLAGS)

ppanalysis: ppanalysis.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(CGAL_LDFLAGS)

ppanalysis.o: ppanalysis.cpp
	$(CXX) $(CXXFLAGS) $(CGAL_CXXFLAGS) -c -o $@ $<

sersic: sersic.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
