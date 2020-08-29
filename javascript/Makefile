PAPAYA2_SRC = ..
CXXFLAGS += -I $(PAPAYA2_SRC)
all:
	emcc -O2 $(CXXFLAGS) --bind -o jspaya.js jspaya.cpp -s WASM=0 -s ALLOW_MEMORY_GROWTH=1 --memory-init-file 0
