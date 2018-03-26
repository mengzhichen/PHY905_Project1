
CXX := c++

CPPFLAGS := -std=c++14 -O3 -larmadillo -lm -llapack -lblas


all: sun_earth.x SEJ.x
	@echo "Compilation complete."

sun_earth.x: Sun_earth.cpp
	$(CXX) -o sun_earth.x Sun_earth.cpp $(CPPFLAGS)

SEJ.x: Sun_earth_jupiter.cpp 
	$(CXX) -o SEJ.x Sun_earth_jupiter.cpp $(CPPFLAGS)


clean:
	-rm *.x *.o
