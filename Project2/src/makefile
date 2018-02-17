
CXX := clang++

CPPFLAGS := -std=c++14 -O3 -larmadillo -lm -llapack -lblas


all: arma cyclic classic
	@echo "Compilation complete."

arma:
	$(CXX) -o armaDiag.x ArmaDiag.cpp $(CPPFLAGS)

cyclic:
	$(CXX) -o cyclicjacobi.x CyclicJacobiDiag.cpp $(CPPFLAGS)

classic:
	$(CXX) -o classicjacobi.x ClassicJacobiDiag.cpp $(CPPFLAGS)

clean:
	-rm *.x *.o
