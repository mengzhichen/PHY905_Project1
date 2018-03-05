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

// This function initializes A just for writing to file 
void Initialization(mat &A, int n){
    for (int i = 0; i < n - 1; i++){
        for (int j = 0; j < n - 1; j++){
            if (i == j)  A(i,j) = 2.0;
            else if (abs(i - j) == 1) A(i,j) = -1;
        }
    } 
}

// This function calculats eigenpolynomials up to rank n
vec eigpoly(vec a, vec b, int n, double x){ 
    vec eig(n);
    vec use(2);
    int sign = 0.0;
    eig(0) = 1.0;
    eig(1) = a(0) - x;
    if (eig(1) <= 0 ) sign = 1;
    if (n>2){
        for (int i = 2; i < n;i++){
             eig(i) = (a(i-1) - x)*eig(i-1) - b(i-2)*b(i-2)*eig(i-2);
            if(eig(i)*eig(i-1)<0 || eig(i)==0) sign+=1;
        }
    }
    use(0) = sign;
    use(1) = eig(n-1);
    return use;
}

// This function uses bisection method to determine subintervals iteratively
void get_subinterval(vec a, vec b, int n, double rmin, double rmax, vec pmin, vec pmax, vec &subint){
    double rmid1 = (rmin+rmax)/2;
    double rmid2 = rmid1;
    vec pmid1 = eigpoly(a, b, n, rmid1);
    vec pmid2 = pmid1;
    //while (pmid1(0)-pmin(0)>1){
    if (pmid1(0)-pmin(0)>1 || pmax(0)-pmid2(0) == 0){
       get_subinterval(a, b, n, rmin, rmid1, pmin, pmid1, subint);
    }
    if (pmax(0)-pmid2(0)>1 || pmid1(0)-pmin(0)==0){
    //while (pmax(0)-pmid2(0)>1){
        get_subinterval(a, b, n, rmid2, rmax, pmid2, pmax, subint);
    }
    if (pmid1(0)-pmin(0)==1){
        subint(pmid1(0)) = rmid1;
        subint(pmin(0)) = rmin;
    }
    if (pmax(0)-pmid2(0)==1){
        subint(pmax(0)) = rmax;
        subint(pmid2(0)) = rmid2;
    }
}

// This function uses bisection method to determine each eigenvalues
vec get_eigenvals(vec a, vec b, int n){
    vec subint(n);
    vec eigval(n-1);
    double tol = 0.00001;
    if (n == 1) eigval(n-1) = a(0);
    double rmin = 0.0 , rmid = 2.0, rmax = 4.0; 
    subint(0) = rmin;
    subint(n-1) = rmax;
    vec pmid(2);   
    vec pmin(2), pmax(2);
    pmin = eigpoly(a, b, n, rmin); 
    pmax = eigpoly(a, b, n, rmax);
    double diff = pmax(0) - pmin(0);
    get_subinterval(a, b, n, rmin, rmax, pmin, pmax, subint);
    for (int i=0; i < n-1; i++){
        rmin = subint(i);
        pmin = eigpoly(a, b, n, rmin);
        if (pmin(0)==0) rmin +=0.000001;
        pmin = eigpoly(a, b, n, rmin);
        rmax = subint(i+1);
        pmax = eigpoly(a, b, n, rmax);
        if (pmax(0)==0) rmax -=0.000001;
        pmax = eigpoly(a, b, n, rmax);
        while (rmax - rmin > tol){
            rmid = (rmin + rmax)/2;
            pmid = eigpoly(a, b, n, rmid);
            if (pmin(1) * pmid(1) < 0) {
                rmax = rmid;
                pmax = pmid;
             }
           else if (pmid(1) == 0){
                break;
             }
            else {
                rmin = rmid;
                pmin = pmid;
            }
            eigval(i) = rmid;
        }
        eigval(i) = rmid;
    }
    return eigval;
}

// This function find the eigenvector for each eigenvalue by solving linear equations
mat get_eigenvecs(vec a, vec b, vec u, int n){
    mat v(n-1,n-1);
    vec a0(n-1);
    double temp=0.0, norm=0.0;
    for (int i = 0; i < n-1; i++){
        for (int k = 0; k < n-1; k++){
            a0(k) = a(k) - u(i);
        }
        norm=0.0;
        v(i,0) = 1.0;
        v(i,1) = -a0(0)/b(0);
        norm = v(i,0) * v(i,0) + v(i,1) * v(i,1);
        for (int j = 2; j < n-1; j++){
            temp = b(j-1);
            v(i,j) = -(b(j-2)*v(i,j-2)/temp + a0(j-1)*v(i,j-1)/temp);
            norm += v(i,j)*v(i,j);
        }
        for (int j = 0; j < n-1; j++){
            v(i,j) = v(i,j)/sqrt(norm);
        }
    }
    return v;
}

////////////////// Main function diagonalize tridiagonal Toeplitz matrix by using the
////////////////// bisection method and compare with analytic eigenvalues.
int main(int argc, char** argv)
{
    int exponent = 0;
    if (argc < 2){
        cout<<"More argument needed"<<endl;
        return 1;
    }
    string filename = "Bisection";
    double pi = acos(-1.0);
    for(int i = 1; i < argc; i++){
        int n = atoi(argv[i]) + 1;             
        double *an_ev = new double[n-1];               //  Vector for analytic egienvalues
        mat eigvec(n-1,n-1);                           //  eigenvectors matrix
        vec a(n-1), b(n-2);                            //  Input Toeplitz matrix
        vec eigval(n-1);  
        mat A(n-1,n-1);              
        
        A.zeros();
        Initialization(A, n);
        for (int i = 0; i < n-2; i++){
            a(i) = 2;
            b(i) = -1;
        }
        a(n-2) = 2;

        //Bisection Diagonalization
        clock_t start = clock();
        eigval = get_eigenvals(a, b, n);
        eigvec = get_eigenvecs(a, b, eigval, n);
        clock_t ends = clock();
        double elapsed_time = (double)(ends - start)/ CLOCKS_PER_SEC * 1000.0;
        cout<<"Execution time: "<<elapsed_time<<endl;
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
        //Write results into the file
        Outputfile(filename, A, eigvec, eigval, n-1, elapsed_time);
        delete [] an_ev;
    }
    return 0;
}
