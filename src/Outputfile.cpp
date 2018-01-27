#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;
ofstream ofile;

inline double u_exact(double x){return 1-(1-exp(-10))*x- exp(-10*x);}

void Outputfile(string filename, int n, double time, double *u, double *x){
    double *rel_err = new double[n+1];         //array of relative error
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

