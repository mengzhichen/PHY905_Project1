
CXX := clang++

CPPFLAGS := -std=c++14 -O3 -larmadillo -lm -llapack -lblas


all: arma cyclic classic one
	@echo "Compilation complete."

arma:
	$(CXX) -o armaDiag.x ArmaDiag.cpp $(CPPFLAGS)

cyclic:
	$(CXX) -o cyclicjacobi.x CyclicJacobi.cpp $(CPPFLAGS)

classic:
	$(CXX) -o classicjacobi.x ClassicJacobi.cpp $(CPPFLAGS)

one:
	$(CXX) -o oneelectron.x ClassicJacobi_one_electron.cpp $(CPPFLAGS)

clean:
	-rm *.x *.
