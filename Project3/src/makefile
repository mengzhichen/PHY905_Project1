
CXX := c++

CPPFLAGS := -std=c++14 -O3 -larmadillo -lm -llapack -lblas


all: Sunearth
	@echo "Compilation complete."

Sunearth:
	$(CXX) -o sun_earth.x Sun_earth.cpp $(CPPFLAGS)

clean:
	-rm *.x *.o
