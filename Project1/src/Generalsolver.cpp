#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include "Outputfile.hpp"

using namespace std;

inline double f(double x){return 100.0*exp(-10*x);} 

int main (int argc, char* argv[])
{
    int exponent = 0;
    if (argc < 2){
        cout<<"More argument needed";
        return 1;
    }
    string filename = "Generalsolution";
    for(int i = 1; i < argc; i++){
        int n = atoi(argv[i]);               //Number of discrete points
        double h = 1.0/n;
        // Setup vectors
        double *a = new double[n+1];         //diagonal elements
        double *b_up = new double[n+1];      //nondiagonal symmetric elements
        double *b_down = new double[n+1];    //nondiagonal symmetric elements
        double *y = new double[n+1];         //y vector in 'Au=y'
        double *u = new double[n+1];         //u vector in 'Au=y'
        double *x = new double[n+1];         //vector of discretized x axis
        //Initialization vectors values
        u[0] = u[n] = 0;
        for(int i = 1; i < n; i++){
            a[i] = 2;
            b_up[i] = b_down[i] = -1;
            x[i] = i*h;
            y[i] = h*h*f(i*h);
        }
        // Forward Substitution
        clock_t start = clock();
        for(int i = 2; i < n; i++){
            double factor = b_down[i]/a[i-1];
            a[i] -= factor * b_up[i-1];
            y[i] -= factor * y[i-1];
        }
        // Bacward Substitution
        u[n-1] = y[n-1]/a[n-1];
        for(int i = n - 2; i > 0; i--){
            u[i] = (y[i] - b_up[i]*u[i+1])/a[i];
        }
        clock_t ends = clock();
        double elapsed_time = (double)(ends - start)/ CLOCKS_PER_SEC * 1000.0;
        
        Outputfile("General_solution", n, elapsed_time, u, x);
        delete [] a; delete [] b_up; delete [] b_down; delete [] y; delete [] u; delete [] x;
    }
    return 0;
}
