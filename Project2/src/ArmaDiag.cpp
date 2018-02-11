#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;
ofstream ofile;

/////////////////// Function for result output.
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

////////////////// Main function diagonalize tridiagonal Toeplitz matrix by using function
////////////////// in armadillo and compare with analytic eigenvalues.
int main(int argc, char** argv)
{
    int exponent = 0;
    if (argc < 2){
        cout<<"More argument needed";
        return 1;
    }
    string filename = "Armadiag";
    double pi = 3.14159265359;
    for(int i = 1; i < argc; i++){
        int n = atoi(argv[i]);                         //Number of discrete points
        double *an_ev = new double[n-1];               //Vector for analytic egienvalues
        mat A(n-1,n-1), eigvec;                        //Input Toeplitz and eigenvectors matrix
        vec eigval;                                    //Vector for armadillo eigenvalues
        
        //Initialization
        A.zeros();
        for (int i = 0; i < n - 1; i++){
             for (int j = 0; j < n - 1; j++){
                 if (i == j) A(i,j) = 2.0;
                 else if (abs(i - j) == 1) A(i,j) = -1;
             }
        }
        
        //Analytical solutions
        cout<<"Analytic eigenvalues:"<<endl;
        for (int i = 0; i < n-1; i++){
            an_ev[i] = 2 * (1 - cos((i+1)*pi/n));
            cout << an_ev[i] <<endl;
        }
        
        //Diagonalization using armadillo lib; time collected.
        clock_t start = clock();
        eig_sym(eigval, eigvec, A);
        clock_t ends = clock();
        double elapsed_time = (double)(ends - start)/ CLOCKS_PER_SEC * 1000.0;
        
        //Results on screen
        A.print("A:");
        eigvec.print("eigvec:");
        eigval.print("eigval:");
        
        //Output results into file
        Outputfile(filename, A, eigvec, eigval, n-1, elapsed_time);
        
        delete [] an_ev;
    }
    return 0;
}
