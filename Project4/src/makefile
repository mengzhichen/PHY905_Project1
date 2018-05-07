
CXX := c++

CPPFLAGS := -std=c++14 -O3 -larmadillo -lm -llapack -lblas


all: vmc.x vmc2.x
	@echo "Compilation complete."

vmc.x: vmc.cpp
	$(CXX) -o vmc.x vmc.cpp $(CPPFLAGS)

vmc2.x: vmc2.cpp
	$(CXX) -o vmc2.x vmc2.cpp $(CPPFLAGS)

clean:
	-rm *.x *.o
