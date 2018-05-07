#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include <random>
#include <functional>
#include "armadillo"

using namespace std;
using namespace arma;

ofstream ofile;



// This function write results into file.
void Outputfile(string filename, double alpha, double beta, double omega, mat LE, mat LE2){
    vec alpha2 = zeros<vec>(20), beta2 = zeros<vec>(20);
    alpha2(0) = alpha+0.025;
    beta2(0) = beta+0.01;
    for (int i=1; i < 20; i++){
        alpha2(i) = alpha2(i-1)+0.025;
        beta2(i) = beta2(i-1)+0.01;
    }
    string fileout = filename;
    fileout.append(".txt");
    ofstream file(fileout);
    file << "The frequency omega = "<< setw(2) << omega << endl;
    cout << setiosflags(ios::fixed);
    file << "The variation parameter alpha = "<< endl;
    alpha2.save(file,raw_ascii);
    file << "The variation parameter beta = "<< endl;
    beta2.save(file,raw_ascii);
    file << "The lowest energies are:" << endl;
    LE.save(file,raw_ascii);
    file << "The variances are: " << endl;
    LE2.save(file,raw_ascii);
    file.close();
}


// This function calculate wave-functions
double Wave_function(double alpha, double beta, double omega, vec r1, vec r2){
    double wf = 0, rr = 0, dr2 = 0, dr = 0, jas = 0;
    for (int i = 0; i < 3; i++){
        rr += r1(i)*r1(i) + r2(i)*r2(i);
        dr2 += (r1(i) - r2(i)) * (r1(i) - r2(i));
    }
    dr = sqrt(dr2);
    jas = dr/(1.0+beta*dr)/2.0;
    wf = exp(-0.5*alpha*omega*rr + jas);
    return wf;
}

// This function calculate the local energy
double Local_energy(double alpha, double beta, double omega, vec r1, vec r2){
    double le = 0, rr = 0, dr2 = 0, dr = 0, jas = 0, jas2 = 0;
    for (int i = 0; i < 3; i++){
        rr += r1(i)*r1(i) + r2(i)*r2(i);
        dr2 += (r1(i) - r2(i)) * (r1(i) - r2(i));
    }
    dr = sqrt(dr2);
    jas = 1.0/(1.0+beta*dr)/2.0;
    jas2 = 2*jas*jas;
    le = 0.5*omega*omega*rr*(1-alpha*alpha) + 3.0*alpha*omega + 1.0/dr + jas2*(alpha*omega*dr - jas2 - 2.0/dr + 4.0*beta*jas) ;
    return le;
}

// This function performs Monte-Carlo sampling

void Monte_carlo(double alpha, double beta, double omega, vec r1, vec r2, mat &LE, mat &LE2, mat &R12){
    int  var = 0, iter = 0, account = 0, account2 = 0, max_var = 20, max_iter = 10000000;
    double wfold, wfnew, le, le2, le_temp, shift, step = 1, per = 0, rr = 0, dr=0, dr_temp=0;
    srand((unsigned)time(NULL));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    
    vec r1new = r1, r2new = r2;
    for (var = 0; var < max_var; var++){
        alpha +=0.025;
        beta = 0.1;
        cout << "beta = " << beta << endl;
        for (int j = 0; j < 20; j++){
            beta += 0.01;
            dr = le = le2 = le_temp = 0;
            for (int i = 0; i < 3; i++){
                r1(i) = step * (dis(gen) - 0.5);
                r2(i) = step * (dis(gen) - 0.5);
            }
            cout << "beta = " << beta << endl;
            wfold = Wave_function(alpha, beta, omega, r1, r2);
        
            account = 0;
            account2 = 0;
            for (iter = 0; iter < max_iter; iter++){
                for (int i = 0; i < 3; i++){
                    r1new(i) = r1(i) + step * (dis(gen) - 0.5);
                    r2new(i) = r2(i) + step * (dis(gen) - 0.5);
                }
                wfnew = Wave_function(alpha, beta, omega, r1new, r2new);
            // Metropolis
                if (dis(gen) <= wfnew*wfnew/(wfold*wfold) ){
                    for (int i = 0; i < 3; i++){
                        r1(i) = r1new(i);
                        r2(i) = r2new(i);
                    }
                    dr_temp = sqrt((r1(0)-r2(0))*(r1(0)-r2(0))+(r1(1)-r2(1))*(r1(1)-r2(1))+(r1(2)-r2(2))*(r1(2)-r2(2)));
                    wfold = wfnew;
                    account +=1;
                }
                if (iter >= 99*max_iter/100){
                    le_temp = Local_energy(alpha, beta, omega, r1, r2);
                    le += le_temp;
                    le2 += le_temp * le_temp;
                    dr += dr_temp;
                    account2+=1;
                }
            }
            per = (double)account/(double)max_iter;
            cout << "alpha = " << alpha << "account = " << account << endl;
            cout << "percentage = " << per << endl;
            cout << "account2 = " << account2 << endl;
            R12(var,j) = dr/account2;
            R12.print("R12:");
            LE(var,j) = le/account2;
            LE.print("LE:");
            LE2(var,j) = le2/account2 - LE(var,j)*LE(var,j);
            LE2.print("LE2:");
        }
    }
}



int main(int argc, char** argv)
{
    if (argc < 2){
        cout<<"More argument needed" <<endl;
        return 1;
    }
    double omega = atof(argv[1]);
    cout<< "omega = " << omega<<endl;
//    double alpha = 0.65, beta = 0.1;
    double alpha = 1.025, beta = 0.18;
    string filename = "Jastrow";
    // 3-D coordinate
    vec r1 = zeros<vec>(3), r2 = zeros<vec>(3);
    mat LE = zeros<mat>(20,20), LE2 = zeros<mat>(20,20), R12 = zeros<mat>(20,20);
    Monte_carlo(alpha, beta, omega, r1, r2, LE, LE2, R12);
    //Outputfile(filename, alpha, beta, omega, LE, LE2);
    return 0;
}
