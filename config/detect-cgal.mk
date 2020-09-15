PRODUCTS += .cgal.mk
MAKEFILES += .cgal.mk

CGD_TEST_COMPILE = $(CXX) $(CXXFLAGS) -o/dev/null 2>/dev/null -c
DUMMY_VARIABLE := $(shell [ -r .cgal.mk ] || \
    ($(CGD_TEST_COMPILE) -std=c++14 ../config/test-recent-cgal.cpp && echo Detected compatible CGAL version: 5.0 and above. >&2 && cp ../config/recent-cgal.mk .cgal.mk) || \
    ($(CGD_TEST_COMPILE) ../config/test-any-cgal.cpp && echo Detected compatible CGAL version: prior to 5.0, we have to link to it. >&2 && cp ../config/older-cgal.mk .cgal.mk) || \
    (echo Found no CGAL headers. >&2 && cp ../config/no-cgal.mk .cgal.mk) \
)
include .cgal.mk
