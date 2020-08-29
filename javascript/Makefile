PAPAYA2_SRC = ../include
CXXFLAGS += -I $(PAPAYA2_SRC)
all:
	emcc -O2 $(CXXFLAGS) --bind -o papaya2.js jspaya2.cpp -s WASM=0 -s ALLOW_MEMORY_GROWTH=1 --memory-init-file 0

clean:
	rm -f papaya2.js
