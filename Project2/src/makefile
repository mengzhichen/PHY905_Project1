
CXX := c++
CPPFLAGS := -std=c++14 -O3 -larmadillo -lm -llapack -lblas

all: arma cyclic classic HO one
	@echo "Compilation complete."

arma: 
	$(CXX) -o armadiag.x ArmaDiag.cpp $(CPPFLAGS)

cyclic:
	$(CXX) -o cyclicjacobi.x CyclicJacobi.cpp $(CPPFLAGS)

classic:
	$(CXX) -o classicjacobi.x ClassicJacobi.cpp $(CPPFLAGS)

HO:
	$(CXX) -o HO.x HO.cpp $(CPPFLAGS)

one:
	$(CXX) -o oneelectron.x ClassicJacobi_one_electron.cpp $(CPPFLAGS)

clean:
	-rm *.x *.
