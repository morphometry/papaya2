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

PRODUCTS = \
    runtests \

default: all
	./runtests

include ../config/shared.mk

all: .ts.mk.hpp runtests

runtests: $(TEST_OBJECTS) runtests.o
	$(CXX) $(CXXFLAGS) -o $@ $(TEST_OBJECTS) runtests.o $(LDFLAGS)
