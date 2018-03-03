#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include "armadillo"

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
        for (int j = i+1; j < n-1; j++){
            norm += A(i,j) * A(i,j);
        }
    }
}

// This function finds largest absoulte element and off-diagobal Frobenius norm
// of input symmetric matrix A.
void maxoffele(mat A, int &r, int &c, double &offnorm, int n){
    offnorm = 0.0;
    double max = 0.0;
    for (int i = 0; i < n-1; i++){
        for (int j = i+1; j < n-1; j++){
            double max_temp = A(i,j) * A(i,j);
            offnorm += max_temp;
            if (max_temp > max){
                max = max_temp;
                r = i;
                c = j;
            }
        }
    }
    offnorm = 2.0 * offnorm;
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
    
// Mainfunction diagonalzie tridiagonal Toeplitz matrix with size n-1 by using classic jacobi iteration.
int main(int argc, char** argv)
{
    int exponent = 0;
    if (argc < 2){
        cout<<"More argument needed";
        return 1;
    }
    string filename = "ClassicJacobi";
    
    for(int i = 1; i < argc; i++){
        int n = atoi(argv[i])+1, iter = 0, maxiter = 1000;
        int r = 0, c = 0;                 // cos(theta) and sin(theta)
        double tol = 0.001;               // error tolerence
         // normal and off-diagonal Frobenius norm of input matrix A
        double norm = 0.0, offnorm = 10.0;
        
        // Define matrices: A and A0 for tridiagonal Toeplitz matrix A; V for eigenvectors.
        mat A0 = zeros<mat>(n-1,n-1);
        mat A = zeros<mat>(n-1,n-1);
        mat V = eye<mat>(n-1,n-1);
        vec eigval(n-1);
        //Initialization
        A.zeros();
        Initialization(A, norm, n);
        Initialization(A0, norm, n);
        // Define global error tolerence
        double delta = tol * norm;
    
        clock_t start = clock();
        
        // The Loop stops until smaller than tolerence error or reach max iter times.
        while (offnorm > delta && iter < maxiter){
            maxoffele(A, r, c, offnorm, n);
            jacobi(A, V, r, c, n);
            iter += 1;
        }
        cout<< "iter = " << iter << endl;
        clock_t ends = clock();
        double elapsed_time = (double)(ends - start)/ CLOCKS_PER_SEC * 1000.0;
        
        ///// place eigenvalues and eigenvecotrs in sequence
        for (int r = 0; r < n - 1; r++){
            for (int c = r; c < n - 1; c++){
                if (A(r,r) > A(c,c)){
                    double A_temp = A(c,c);
                    A(c,c) = A(r,r);
                    A(r,r) = A_temp;
                    V.swap_cols(r,c);
                }
            }
        }
        for (int i = 0; i < n-1; i++){
            eigval(i) = A(i,i);
        }
        A.print("A:");
        V.print("eigenvecors:");
        Outputfile(filename, A0, V, eigval, n-1, elapsed_time);
    }
    return 0;
}
