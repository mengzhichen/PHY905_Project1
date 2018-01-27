
CXX := g++

CPPFLAGS := -std=c++11 -O3 
CPPFLAGS2 := -std=c++11 -O3 -larmadillo -lm


all: general special lu
	@echo "Compilation complete."

general:
	$(CXX) -c -o Outputfile.o Outputfile.cpp $(CPPFLAGS)
	$(CXX) -o general.x Generalsolver.cpp Outputfile.o $(CPPFLAGS)

special:
	$(CXX) -c -o Outputfile.o Outputfile.cpp $(CPPFLAGS)
	$(CXX) -o special.x Specialsolver.cpp Outputfile.o $(CPPFLAGS)

lu:
	$(CXX) -c -o Outputfile.o Outputfile.cpp $(CPPFLAGS)
	$(CXX) -o lu.x LUdecomposition.cpp Outputfile.o $(CPPFLAGS2)


clean:
	-rm *.x
