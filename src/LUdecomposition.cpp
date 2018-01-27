#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;
ofstream ofile;

inline double f(double x){return 100.0*exp(-10*x);}
inline double u_exact(double x){return 1-(1-exp(-10))*x- exp(-10*x);}

void Outputfile(int n, double time, double *u, double *x){
    cout << n <<endl;
    double *rel_err = new double[n+1];         //array of relative errors
    string filename = "LU_soulution";
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

int main(int argc, char** argv)
{
    int exponent = 0;
    if (argc < 2){
        cout<<"More argument needed";
        return 1;
    }
    string filename = "LU_decomposition";
    
    for(int i = 1; i < argc; i++){
        cout<< i<<endl;
        int n = atoi(argv[i]);                        //Number of discrete points
        double h = 1.0/n;
        double *y = new double[n];                  //y vector in 'Au=y'
        double *u = new double[n];                  //u vector in 'Au=y'
        double *y_temp = new double[n-1];             //vector satisfies L*y_temp=y
        double *x = new double[n];                  //vector of discretized x axis
        mat A(n-1,n-1), Low(n-1,n-1), Up(n-1,n-1);    //Define matrices
        //Initialization
        A.zeros();
        for (int i = 0; i < n - 1; i++){
             for (int j = 0; j < n - 1; j++){
                 if(i == j) A(i,j) = 2.0;
                 else if (abs(i - j) == 1) A(i,j) = -1;
             }
            x[i] = i*h;
            y[i] = h*h*f(i*h);
        }
        x[n-1] = (n-1)*h;
        y[n-1] = h*h*f((n-1)*h);
        
        clock_t start = clock();
        //LU decomposition: A = L*U
        lu(Low,Up,A);
        //Forward Substitution
        y_temp[0] = y[1];
        for (int i = 1; i < n - 1; i++){
            y_temp[i] = y[i+1];
            for (int j = 0; j < i; j++){
                y_temp[i] -= Low(i,j)*y_temp[j];
            }
        }
        //Backward Substitution
        u[n-1] = y_temp[n-2]/Up(n-2,n-2);
        for (int i = n-2; i > 0; i--){
            u[i] = y_temp[i-1];
            for (int j = i; j <= n-2; j++){
                u[i] -= Up(i-1,j)*u[j+1];
            }
            u[i] = u[i]/Up(i-1,i-1);
        }
        clock_t ends = clock();
        double elapsed_time = (double)(ends - start)/ CLOCKS_PER_SEC * 1000.0;
        
        Outputfile(n, elapsed_time, u, x);
        delete [] u; delete [] y; delete [] y_temp; delete [] x;
    }
    return 0;
}
