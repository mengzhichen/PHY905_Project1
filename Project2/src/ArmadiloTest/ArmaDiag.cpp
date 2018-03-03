#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include "armadillo"
#include "ArmaDiag.h"

using namespace std;
using namespace arma;

////////////////// Main function diagonalize tridiagonal Toeplitz matrix by using function
////////////////// in armadillo.
vec Toeplitz(int n)
{
        mat A(n,n), eigvec;                        //Input Toeplitz and eigenvectors matrix
        vec eigval;                                    //Vector for armadillo eigenvalues
        
        //Initialization
        A.zeros();
        for (int i = 0; i < n; i++){
             for (int j = 0; j < n; j++){
                 if (i == j) A(i,j) = 2.0;
                 else if (abs(i - j) == 1) A(i,j) = -1;
             }
        }
        
        //Diagonalization using armadillo lib; time collected.
        clock_t start = clock();
        eig_sym(eigval, eigvec, A);
        clock_t ends = clock();
        double elapsed_time = (double)(ends - start)/ CLOCKS_PER_SEC * 1000.0;

    return eigval;
}

vec Analytical(int n)
{
    double pi = 3.14159265359;
    vec ana_eigval(n);
    for (int i = 0; i < n; i++){
        ana_eigval(i) = 2 * (1 - cos((i+1)*pi/(n+1)));
    }
    return ana_eigval;
}
