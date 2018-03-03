#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include "armadillo"
#include "ClassicJacobi.h"

using namespace std;
using namespace arma;

// This function write results into file.
void Outputfile(string filename, mat A, mat B, vec v, int n, double time){
    string fileout = filename;
    string size = to_string(n);
    fileout.append("_"+size+".txt");
    ofstream file(fileout);
    file << "Diagonalization takes time: " <<time << "ms"  << endl;
    file << "The input matrix size is: "<<size<<"*"<<size<<endl;
    file << "The input Matrix A is:" << endl;
    A.save(file,raw_ascii);
    file << "The eigenvectors are: " << endl;
    B.save(file,raw_ascii);
    file << "The eigenvalues are: " << endl;
    v.save(file,csv_ascii);
    file.close();
}

// This function initializes the matrix A which we are going to diagonalize
void Initialization(mat &A, double &norm, int n){
    for (int i = 0; i < n - 1; i++){
        for (int j = 0; j < n - 1; j++){
            if (i == j)  A(i,j) = 2.0;
            else if (abs(i - j) == 1) A(i,j) = -1;
            norm += A(i,j) * A(i,j);
        }
    }
    
}

// This function calculates Frobenius norm to check it's conservativeness.
// It should only used for test!
void Frobeniusnorm(mat A, double &norm, int n){
    norm = 0.0;
    for (int i = 0; i < n-1; i++){
        for (int j = 0; j < n-1; j++){
            norm += A(i,j) * A(i,j);
        }
    }
}

// This function finds largest absoulte element of input symmetric matrix A.
void maxoffele(mat A, int &r, int &c, double &offmax, int n){
    offmax = 0;
    for (int i = 0; i < n-1; i++){
        for (int j = i+1; j < n-1; j++){
            double max_temp = A(i,j) * A(i,j);
            if (max_temp > offmax){
                offmax = max_temp;
                r = i;
                c = j;
            }
        }
    }
}

// This function performs classical Jacobi's method to remove off-diagonal elements.
void jacobi(mat &A, mat &R, int r, int c, int n){
    double si, co;
    if (A(r,c) != 0){
        double tau, tan;
        tau = (A(c,c) - A(r,r))/2.0/A(r,c);
        if (tau >= 0) {
            tan = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else tan = -1.0/(-tau + sqrt(1.0 + tau*tau));
        co = 1/sqrt(1 + tan *tan);
        si = co * tan;
    }
    else{
        co = 1.0, si = 0.0;
    }
    double A_rr = A(r,r);
    double A_cc = A(c,c);
    A(r,r) = A_rr * co * co - 2.0 * A(r,c) * si * co + A_cc * si *si;
    A(c,c) = A_cc * co * co + 2.0 * A(r,c) * si * co + A_rr * si *si;
    A(r,c) = A(c,r) = 0.0;
    for (int i = 0; i < n-1; i++){
        if (i != r && i != c){
            double A_ir = A(i,r);
            double A_ic = A(i,c);
            A(r,i) = A(i,r) = A_ir*co - A_ic*si;
            A(c,i) = A(i,c) = A_ic*co + A_ir*si;
        }
    /// Calculate Eigenvectors R
    double R_ir = R(i,r);
    double R_ic = R(i,c);
    R(i,r) = R_ir * co - R_ic * si;
    R(i,c) = R_ic * co + R_ir * si;
    }
}
