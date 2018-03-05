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


// This function performs classical Jacobi's method to remove off-diagonal elements
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
//// This function calculates cos(theta) and sin(theta) which are elements in S and S_t.
//void Cyclic_jacobi(mat A, int r, int c, double &co, double &si){
        //double tau, tan;
        //if (A(r,c) != 0){
            //tau = (A(c,c) - A(r,r))/2/A(r,c);
            //if (tau >= 0) {
                //tan = 1.0/(tau + sqrt(1.0 + tau*tau));
            //}
            //else tan = -1.0/(-tau + sqrt(1.0 + tau*tau));
            //co = 1/(sqrt(1 + tan*tan));
            //si = tan * co;
        //}
        //else{
            //co = 1.0, si = 0.0;
        //}
    //}
                            
// Mainfunction diagonalzie tridiagonal Toeplitz matrix with size n-1 by using cyclic jacobi iteration.
int main(int argc, char** argv)
{
    int exponent = 0;
    if (argc < 2){
        cout<<"More argument needed";
        return 1;
    }
    string filename = "CyclicJacobi";
    
    for(int i = 1; i < argc; i++){
        int n = atoi(argv[i]) + 1, iter = 1, maxiter = 100;
        double co = 0.0;
        double si = 0.0;
        double tol = 0.00000001;
        double pi = acos(-1.0);
        double *an_ev = new double[n-1];               //  Vector for analytic egienvalues
        
        // Define matrices: A and A0 for tridiagonal Toeplitz matrix A; V for eigenvectors;
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
        
        //  Calculate diagonal and off-diagonal Frobenius norm
        double norm = 0.0, off_norm = 0.0;
        
        for (int i = 0; i < n - 1; i++){
            for (int j = 0; j < n - 1; j++){
                norm += A(i,j)*A(i,j);
                if (i != j) off_norm += A(i,j)*A(i,j);
            }
        }
        // Define error global tolerence
        double delta = tol * norm;
        
        clock_t start = clock();
        // The loop stops until smaller than tolerence error or reach max iter times
        while (off_norm > delta && iter < maxiter){
            for (int r = 0; r < n - 2; r++){
                for (int c = r + 1; c < n - 1; c++){
                    jacobi(A, V, r, c, n);
                }
            }
            off_norm = norm;
            for (int i = 0; i < n - 1; i++){
                off_norm -= A(i,i)*A(i,i);
            }
            iter += 1;
        }
        cout<<"iter = "<<iter<<endl;
        clock_t ends = clock();
        double elapsed_time = (double)(ends - start)/ CLOCKS_PER_SEC * 1000.0;
        
        // place eigenvalues and eigenvecotrs in sequence
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

        //Analytical solutions
        cout<<"First five Analytic eigenvalues:"<<endl;
        for (int i = 0; i < n-1; i++){
            an_ev[i] = 2 * (1 - cos((i+1)*pi/n));
        }
        cout <<setiosflags(ios::fixed);
        for (int i = 0; i < 5; i++){
            cout << setprecision(5)<<an_ev[i] <<endl;
        }

        cout.precision(5);
        cout.setf(ios::fixed);
        eigval.raw_print(cout,"eigval:");
        cout<<"Execution time: "<<elapsed_time<<endl;
       // A.print("A:");
        //V.print("V:");
        Outputfile(filename, A0, V, eigval, n-1, elapsed_time);
    }
    return 0;
}
