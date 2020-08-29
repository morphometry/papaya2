# defaults
CXXFLAGS += -g -O3
CXXFLAGS += -Wall -Werror
FITS_LDFLAGS = -lCCfits -lcfitsio
CXXFLAGS += -I../include

include ../cgal.mk

# modify defaults here if necessary
-include features.mk

MAKEFILES = \
    Makefile \

BINARIES = \
    banana \
    sersic \
    imganalysis \
    ppanalysis \

default: .ts.mk.hpp imganalysis ppanalysis

all: default banana sersic

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

# special rule, to be able to load FITS containers
banana: banana.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(FITS_LDFLAGS)

ppanalysis: ppanalysis.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(CGAL_LDFLAGS)

ppanalysis.o: ppanalysis.cpp
	$(CXX) $(CXXFLAGS) $(CGAL_CXXFLAGS) -c -o $@ $<

sersic: sersic.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(BINARIES) *.o .ts.mk.hpp

.PHONY: default all clean depend
