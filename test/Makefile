# defaults
CXXFLAGS += -g -O3 -std=c++11
CXXFLAGS += -Wall -Werror
CXXFLAGS += -I../include

MAKEFILES = \
    Makefile \

TEST_OBJECTS = \
    test_image.o \
    test_tools.o \
    test_imt.o \

BINARIES = \
    runtests \

default: all
	./runtests

all: .ts.mk.hpp runtests

# hack to make any code recompile if Makefile changes
.ts.mk.hpp: $(MAKEFILES)
	@touch $@

# automagic dependencies
depend: *.cpp *.hpp .ts.mk.hpp
	$(CXX) $(CXXFLAGS) -MM *.cpp >.depend
-include .depend

runtests: $(TEST_OBJECTS) runtests.o
	$(CXX) $(CXXFLAGS) -o $@ $(TEST_OBJECTS) runtests.o $(LDFLAGS)

clean:
	rm -f $(BINARIES) *.o .ts.mk.hpp

.PHONY: default all clean depend
