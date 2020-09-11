MAKEFILES += \
    Makefile \
    features.mk \

# create features.mk if absent
features.mk:
	@touch $@

clean:
	rm -f $(PRODUCTS) *.o .ts.mk.hpp .depend.mk

# hack to make any code recompile if Makefile changes
.ts.mk.hpp: $(MAKEFILES)
	@touch $@

# automagic dependencies
depend: *.cpp *.hpp .ts.mk.hpp
	$(CXX) $(CXXFLAGS) -MM *.cpp >.depend.mk
-include .depend.mk

.PHONY: default all clean depend
