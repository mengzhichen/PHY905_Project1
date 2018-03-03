CXX= c++
# Here we define the library functions we need
CPPFLAGS = -std=c++14 -O3
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG= testcode.x
${PROG} :	   	tests-main.o  test-functions.o ArmaDiag.o
			${CXX} ${CPPFLAGS} tests-main.o test-functions.o ArmaDiag.o ${LIB} -o ${PROG}

tests-main.o :			tests-main.cpp 
		        	${CXX} ${CPPFLAGS} -c tests-main.cpp

test-functitons.o :			test-functions.cpp
		        	${CXX} ${CPPFLAGS} -c test-functions.cpp

ArmaDiag.o :		ArmaDiag.cpp 
		        	${CXX} ${CPPFLAGS} -c ArmaDiag.cpp

clean:
	-rm *.x *.o
