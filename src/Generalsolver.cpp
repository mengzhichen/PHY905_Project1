#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include <string>

using namespace std;
ofstream ofile;

inline double f(double x){return 100.0*exp(-10*x);}
inline double u_exact(double x){return 1-(1-exp(-10))*x- exp(-10*x);}

void Outputfile(int n, double time, double *u, double *x){
    cout << n <<endl;
    double *rel_err = new double[n+1];         //array of relative errors
    string filename = "Generalsolution";
    string fileout = filename;
    string argument = to_string(n);
    fileout.append("_"+argument+".txt");
    ofile.open(fileout);
    ofile <<setiosflags(ios::showpoint | ios::fixed | ios::right) << setw(20) << "Substitution Time : " <<time << "ms"  << endl;
    cout << "Substitution Time : "<<time <<"ms"<< endl;
    ofile << setw(15) << "x" << setw(15) << "u" << setw(15) << "u_exact" << setw(15) << "rel-error" <<endl;
    double max_rel_err = -100;
    for(int i = 1; i < n; i++){
        double x_val = x[i];
        rel_err[i] = log10(fabs((u_exact(x_val) - u[i])/u_exact(x_val)));
        if (rel_err[i] > max_rel_err){
            max_rel_err = rel_err[i];
        }
        ofile << setw(15) << setprecision(8) << x_val;
        ofile << setw(15) << setprecision(8) << u[i];
        ofile << setw(15) << setprecision(8) << u_exact(x_val);
        ofile << setw(15) << setprecision(8) << rel_err[i] << endl;
    }
    ofile.close();
    cout << "Maxium relative error is: " << max_rel_err << endl;
    delete [] rel_err;
}

int main (int argc, char* argv[])
{
    int exponent = 0;
    if (argc < 2){
        cout<<"More argument needed";
        return 1;
    }
    string filename = "Generalsolution";
    for(int i = 1; i < argc; i++){
        int n = atoi(argv[i]);          //Number of discrete points
        double h = 1.0/n;
        // Setup vectors
        double *a = new double[n+1];      //diagonal elements
        double *b_up = new double[n+1];    //nondiagonal symmetric elements
        double *b_down = new double[n+1];    //nondiagonal symmetric elements
        double *y = new double[n+1];      //y vector in 'Au=y'
        double *u = new double[n+1];      //u vector in 'Au=y'
        double *x = new double[n+1];      //vector of discretized x axis
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
        
        Outputfile(n, elapsed_time, u, x);
        delete [] a; delete [] b_up; delete [] b_down; delete [] y; delete [] u; delete [] x;
    }
    return 0;
}
