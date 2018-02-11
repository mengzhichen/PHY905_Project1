#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include "armadillo"

#define pi 3.14159265359

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

// This function calculates cos(theta) and sin(theta) which are elements in S and S_t.
void Cyclic_jacobi(mat A, int r, int c, double &co, double &si){
        double tau, tan;
        tau = (A(c,c) - A(r,r))/2/A(r,c);
        if (A(r,c) != 0){
            if (tau >= 0){
                tan = -tau + sqrt(1 + tau * tau);}
            else {
                tan = -tau - sqrt(1 + tau * tau);}
            co = 1/(sqrt(1 + tan*tan));
            si = tan * co;
        }
        else{
            co = 1.0, si = 0.0;
        }
    }
                            
// Mainfunction diagonalzie tridiagonal Toeplitz matrix with size n-1 by using cyclic jacobi method.
int main(int argc, char** argv)
{
    int exponent = 0;
    if (argc < 2){
        cout<<"More argument needed";
        return 1;
    }
    string filename = "CyclicJacobi";
    for(int i = 1; i < argc; i++){
        int n = atoi(argv[i]), iter = 1, maxiter = 1000;
        double co = 0.0;
        double si = 0.0;
        double tol = 0.001;
        
        // Define matrices: A for tridiagonal Toeplitz matrix A; V for eigenvectors;
        // S and S_t are orthogonal transoformation matrices.
        mat A0(n-1,n-1), A(n-1,n-1);
        mat V =eye<mat>(n-1,n-1);
        mat S =eye<mat>(n-1,n-1);
        mat S_t =eye<mat>(n-1,n-1);
        vec eigval(n-1);
        //Initialization
        A.zeros();
        for (int i = 0; i < n - 1; i++){
            for (int j = 0; j < n - 1; j++){
                if (i == j) A0(i,j) = A(i,j) = 2.0;
                else if (abs(i - j) == 1) A0(i,j) = A(i,j) = -1;
            }
        }
        
        //Calculate diagonal and off-diagonal Frobenius norm
        double norm = 0.0, off_norm = 0.0;
        
        for (int i = 0; i < n - 1; i++){
            for (int j = 0; j < n - 1; j++){
                norm += A(i,j)*A(i,j);
                if (i != j) off_norm += A(i,j)*A(i,j);
            }
        }
        //Define error global tolerence
        double delta = tol * norm;
        
        //Loop stops until smaller than tolerence error or reach max iter
        clock_t start = clock();
        
        while (off_norm > delta || iter > maxiter){
            for (int r = 0; r < n - 2; r++){
                for (int c = r + 1; c < n - 1; c++){
                    Cyclic_jacobi(A, r, c, co, si);
                    S = S_t = eye<mat>(n-1,n-1);
                    S(r,r) = S(c,c) = S_t(r,r) = S_t(c,c) = co;
                    S(r,c) = S_t(c,r) =  si;
                    S(c,r) = S_t(r,c) = -si;
                    A = S_t*A*S;
                    V = V*S;
                }
            }
            off_norm = 0;
            for (int i = 0; i < n - 1; i++){
                for (int j = 0; j < n - 1; j++){
                    if (i != j) off_norm += A(i,j)*A(i,j);
                }
            }
            cout<<"iter = "<<iter<<endl;
            iter += 1;
        }
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
            eigval(r) = A(r,r);
        }
        A.print("A:");
        V.print("V:");
        Outputfile(filename, A0, V, eigval, n-1, elapsed_time);
    }
    return 0;
}
