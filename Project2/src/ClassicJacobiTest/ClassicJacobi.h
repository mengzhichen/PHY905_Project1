#ifndef JACOBI_H
#define	JACOBI_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

void Outputfile(string, mat, mat, vec, int, double);
void Initialization(mat&, double&, int);
void Frobeniusnorm(mat, double&, int);
void maxoffele(mat, int&, int&, double&, int);
void jacobi(mat&, mat&, int, int, int);

#endif /* JACOBI_H */

