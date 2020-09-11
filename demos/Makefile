# defaults
CXXFLAGS += -g -O3
CXXFLAGS += -Wall -Werror
FITS_LDFLAGS = -lCCfits -lcfitsio
CXXFLAGS += -I../include

# modify defaults here if necessary
-include features.mk

PRODUCTS = \
    banana \
    sersic \
    imganalysis \
    ppanalysis \

default: imganalysis ppanalysis

include ../config/shared.mk
include ../config/detect-cgal.mk

all: default banana sersic

imganalysis: imganalysis.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

# special rule, to be able to load FITS containers
banana: banana.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(FITS_LDFLAGS)

ppanalysis: ppanalysis.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(CGAL_LDFLAGS)

sersic: sersic.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
