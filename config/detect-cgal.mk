PRODUCTS += .cgal.mk
MAKEFILES += .cgal.mk
.cgal.mk:
	@echo Please ignore the above error about .cgal.mk not existing.
	@($(CXX) -o /dev/null -c ../config/test-recent-cgal.cpp 2>/dev/null && echo Detected CGAL 5.0 and above. && cp ../config/recent-cgal.mk .cgal.mk) || \
            ($(CXX) -o /dev/null -c ../config/test-any-cgal.cpp 2>/dev/null && echo Detected CGAL prior to 5.0. && cp ../config/older-cgal.mk .cgal.mk) || \
            (echo Found no CGAL headers. && cp ../config/no-cgal.mk .cgal.mk)
include .cgal.mk
