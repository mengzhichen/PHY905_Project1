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
void Outputfile(string filename, double alpha, double omega, vec LE, vec LE2){
    vec alpha2 = zeros<vec>(20);
    alpha2(0) = alpha+0.025;
    for (int i=1; i < 20; i++){
        alpha2(i) = alpha2(i-1)+0.025;
    }
    string fileout = filename;
    fileout.append(".txt");
    ofstream file(fileout);
    file << "The frequency omega = "<< setw(2) << omega << endl;
    cout << setiosflags(ios::fixed);
    file << "The variation parameter alpha = "<< endl;
    alpha2.save(file,raw_ascii);
    file << "The lowest energies are:" << endl;
    LE.save(file,raw_ascii);
    file << "The variances are: " << endl;
    LE2.save(file,raw_ascii);
    file.close();
}


// This function calculate wave-functions
double Wave_function(double alpha, double omega, vec r1, vec r2){
    double wf = 0, rr = 0;
    for (int i = 0; i < 3; i++){
        rr += r1(i)*r1(i) + r2(i)*r2(i);
    }
    wf = exp(-alpha*omega*rr/2);
    return wf;
}


// This function calculate the local energy
double Local_energy(double alpha, double omega, vec r1, vec r2){
    double le = 0, rr = 0, dr2 = 0, dr = 0;
    for (int i = 0; i < 3; i++){
        rr += r1(i)*r1(i) + r2(i)*r2(i);
        dr2 += (r1(i) - r2(i)) * (r1(i) - r2(i));
    }
    dr = sqrt(dr2);
    le = 0.5*omega*omega*rr*(1-alpha*alpha) + 3.0*alpha*omega + 1/dr;
//    le = 0.5*omega*omega*rr*(1-alpha*alpha) + 3.0*alpha*omega;
    return le;
}

// This function calculate the local energy
double Local_kenergy(double alpha, double omega, vec r1, vec r2){
    double le = 0, rr = 0, dr2 = 0, dr = 0;
    for (int i = 0; i < 3; i++){
        rr += r1(i)*r1(i) + r2(i)*r2(i);
        dr2 += (r1(i) - r2(i)) * (r1(i) - r2(i));
    }
    dr = sqrt(dr2);
    le = -0.5*omega*omega*rr*alpha*alpha + 3.0*alpha*omega;
    return le;
}

// This function calculate the local energy
double Local_penergy(double alpha, double omega, vec r1, vec r2){
    double le = 0, rr = 0, dr2 = 0, dr = 0;
    for (int i = 0; i < 3; i++){
        rr += r1(i)*r1(i) + r2(i)*r2(i);
        dr2 += (r1(i) - r2(i)) * (r1(i) - r2(i));
    }
    dr = sqrt(dr2);
    le = 0.5*omega*omega*rr + 1/dr;
    //le = 0.5*omega*omega*rr;
    return le;
}

// This function performs Monte-Carlo sampling
void Monte_carlo(double alpha, double omega, double step, vec r1, vec r2, vec &LE, vec &LE2, vec &R12, vec &KE, vec &PE){
    int  var = 0, iter = 0, account = 0, account2 = 0, max_var = 20, max_iter = 10000000;
    double wfold, wfnew, le, le2, le_temp, shift, per = 0, rr = 0, dr=0, dr_temp=0;
    double ke_temp, ke, pe_temp, pe;
    srand((unsigned)time(NULL));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    
    vec r1new = r1, r2new = r2;
    for (var = 0; var < max_var; var++){
        alpha +=0.025;
        cout<< "alpha = "<< alpha<<endl;
        //alpha = 1;
        ke=pe=dr = le = le2 = 0;
        for (int i = 0; i < 3; i++){
            r1(i) = step * (dis(gen) - 0.5);
            r2(i) = step * (dis(gen) - 0.5);
        }
        wfold = Wave_function(alpha, omega, r1, r2);
        cout << "wfold = " << wfold << endl;
        
        account = 0;
        account2 = 0;
        for (iter = 0; iter < max_iter; iter++){
            for (int i = 0; i < 3; i++){
                r1new(i) = r1(i) + step * (dis(gen) - 0.5);
                r2new(i) = r2(i) + step * (dis(gen) - 0.5);
            }
            wfnew = Wave_function(alpha, omega, r1new, r2new);
            //cout << "wfnew= " << wfnew << endl;
            // Metropolis
            if (dis(gen) <= wfnew*wfnew/(wfold*wfold) ){
                for (int i = 0; i < 3; i++){
                    r1(i) = r1new(i);
                    r2(i) = r2new(i);
                }
                dr_temp = sqrt((r1(0)-r2(0))*(r1(0)-r2(0))+(r1(1)-r2(1))*(r1(1)-r2(1))+(r1(2)-r2(2))*(r1(2)-r2(2)));
                wfold = wfnew;
                account = account+1;
            }
            if (iter >= 9999*max_iter/10000){
                le_temp = Local_energy(alpha, omega, r1, r2);
                ke_temp = Local_kenergy(alpha, omega, r1, r2);
                pe_temp = Local_penergy(alpha, omega, r1, r2);
                le += le_temp;
                ke += ke_temp;
                pe += pe_temp;
                le2 += le_temp * le_temp;
                dr += dr_temp;
                account2+=1;
            }
        }
        for (int i = 0; i < 3; i++){
            rr += r1(i)*r1(i) + r2(i)*r2(i);
        }
        per = (double)account/(double)max_iter;
        cout.setf(ios::fixed);
        cout << "alpha = " << alpha << "account = " << account << endl;
        cout << "percentage = " << per << endl;
        R12(var) = dr/account2;
        //cout << "R12" << R12 <<endl;
        KE(var) = ke/account2;
        cout <<fixed << setprecision(2)<< "KE" << KE <<endl;
        PE(var) = pe/account2;
        cout  << setprecision(10)<< "PE" << PE <<endl;
        LE(var) = le/account2;
        cout << "LE" << LE <<endl;
        LE2(var) = le2/account2 - LE(var)*LE(var);
        cout << "LE2" << LE2 <<endl;
    }
}



int main(int argc, char** argv)
{
    if (argc < 2){
        cout<<"More argument needed" <<endl;
        return 1;
    }
    double omega = atof(argv[1]);
    double step = atof(argv[2]);
    cout<< "omega = " << omega<<endl;
    double alpha = 0.4;
    string filename = "Simple";
    // 3-D coordinate
    vec r1 = zeros<vec>(3), r2 = zeros<vec>(3), R12 = zeros<vec>(20), LE = zeros<vec>(20), LE2 = zeros<vec>(20), KE = zeros<vec>(20), PE = zeros<vec>(20);
    Monte_carlo(alpha, omega, step, r1, r2, LE, LE2, R12, KE, PE);
    //Outputfile(filename, alpha, omega, LE, LE2);
    return 0;
}
