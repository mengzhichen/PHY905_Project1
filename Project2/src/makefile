
CXX := clang++

CPPFLAGS := -std=c++14 -O3 -larmadillo -lm -llapack -lblas


all: arma jacobi
	@echo "Compilation complete."

arma:
	$(CXX) -o arma.x arma_diag.cpp $(CPPFLAGS)

jacobi:
	$(CXX) -o jacobi.x Jacobi.cpp $(CPPFLAGS)

clean:
	-rm *.x *.o
