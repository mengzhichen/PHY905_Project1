
CXX := c++
CPPFLAGS := -std=c++14 -O3 -larmadillo -lm -llapack -lblas

all: arma cyclic classic oneelectron onejacobi bisection twoelectron
	@echo "Compilation complete."
.PHONY: all

arma: ArmaDiag.cpp
	$(CXX) -o armadiag.x ArmaDiag.cpp $(CPPFLAGS)

cyclic: CyclicJacobi.cpp
	$(CXX) -o cyclicjacobi.x CyclicJacobi.cpp $(CPPFLAGS)

classic: ClassicJacobi.cpp
	$(CXX) -o classicjacobi.x ClassicJacobi.cpp $(CPPFLAGS)

oneelectron: oneelectron_arma.cpp
	$(CXX) -o oneeletron.x oneelectron_arma.cpp $(CPPFLAGS)

onejacobi: ClassicJacobi_one_electron.cpp
	$(CXX) -o oneelectron_jacobi.x ClassicJacobi_one_electron.cpp $(CPPFLAGS)

bisection: Bisection.cpp
	$(CXX) -o bisection.x Bisection.cpp $(CPPFLAGS)

twoelectron: twoelectron_arma.cpp
	$(CXX) -o twoelectron.x twoelectron_arma.cpp $(CPPFLAGS)


clean:
	-rm *.x *.
