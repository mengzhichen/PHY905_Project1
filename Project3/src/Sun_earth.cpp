#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

ofstream ofile;
// function declarations
void Outputfile(string, int, double, double, vec, vec, vec, vec, vec, vec, vec, vec, vec, vec);
void Verlet( int, double, double);
void Euler( int, double, double);

// This function write results into file.
void Outputfile(string filename, int num_steps, double exet, double t_end, vec t, vec x, vec y, vec r, vec vx, vec vy, vec kin, vec pot, vec tot, vec Lz){
    string fileout = filename;
    string size = to_string(t_end);
    fileout.append("_"+size+".txt");
    ofstream file(fileout);
    file << "The excution takes:"<<setw(2) << exet <<setw(1)<<"ms"<<endl;
    file << "time(yr) " << setw(5) << "x" << setw(10) << "y" << setw(10) << "r" << setw(10) << "vx"<< setw(10) << "vy"  << setw(13) << "kin_ene" << setw(10) << "pot_ene"<< setw(10) << "tot_ene" << setw(8) << "Lz"<< endl;
    cout << setiosflags(ios::fixed);
    for (int i = 0; i < num_steps; i++){
        file << fixed << setprecision(6) << t(i) << setw(10) << x(i) << setw(10) << y(i) << setw(10) << r(i) << setw(10) << vx(i) << setw(10) << vy(i) << setw(10) << kin(i) << setw(10) << pot(i) << setw(10) << tot(i) << setw(10) << Lz(i) << endl;
    }
    file.close();
}

// This function calculates time-evolution by velovity Verlet method.
void Verlet(int num_steps, double t_end, double pi){
    double h = t_end /num_steps;                       //time steps
    double fac = 4.0 * h * pi * pi;                    //prefactor for acceleration
    double m_ear = 0.000003;                           //earth mass (in unit of M_sun(Sun mass) = 1)
    double G = 4*pi*pi;                                //Gravatational constant(in unit of AU^3*M_sun^-1*yr^-2 )
    double kin_fac = 0.5*m_ear, pot_fac = -G*m_ear;    //prefactors for kinetic and potential energy
    num_steps += 1;
    //Define and initialize coordinates x, y and their speed vx and vy
    vec x(num_steps), y(num_steps), vx(num_steps), vy(num_steps), t(num_steps), r(num_steps);
    x(0) = 1.0; y(0) = 0.0; vx(0) = 0.0;  vy(0)=2.0*pi; t(0) = 0.0;
    r(0) = 1.0;
    //Define and initialize kinetic, potential and total energy (in unit of M_sun*AU^2*yr^-2),
    //as well as angular momentum Lz (in unit of M_sun*AU^2*yr^-1)
    vec kin_ene(num_steps), pot_ene(num_steps), tot_ene(num_steps), Lz(num_steps);
    kin_ene(0) = kin_fac*vy(0)*vy(0); pot_ene(0) = pot_fac/r(0); tot_ene(0) =  kin_ene(0)+pot_ene(0);
    Lz(0) = m_ear*vy(0);
    //acceleration times h at step i and i+1;
    double axh0 = 0.0, ayh0 = 0.0, axh1 = 0.0, ayh1 = 0.0;
//    double emax = 1.5, emin = 1.3;
//    double emid = 1.4;
    
//    double s = 1;
//    int jj = 0;
    axh0 = fac * x(0)/(r(0)*r(0)*r(0));
    ayh0 = fac * y(0)/(r(0)*r(0)*r(0));
    clock_t start = clock();
//    while (s < num_steps - 5){
//        vy(0)=emid*2.0*pi;
//        jj = 0;
//        cout << "emid = "<<emid<<endl;
        for (int i = 1; i < num_steps; i++){
        //Calculate time-evolution of coordinates and velocities
            x(i) = x(i-1) + h * vx(i-1) - h * axh0/2;
            y(i) = y(i-1) + h * vy(i-1) - h * ayh0/2;
            r(i) = sqrt(x(i) * x(i) + y(i) * y(i));
            axh1 = fac * x(i)/(r(i)*r(i)*r(i));
            ayh1 = fac * y(i)/(r(i)*r(i)*r(i));
            vx(i) = vx(i-1) - 0.5*(axh0+axh1);
            vy(i) = vy(i-1) - 0.5*(ayh0+ayh1);
            axh0 = axh1;
            ayh0 = ayh1;
        //Calculate time-evolution of energies and angular momentum
            kin_ene(i) = kin_fac*(vx(i)*vx(i) + vy(i)*vy(i));
            pot_ene(i) = pot_fac/r(i);
            tot_ene(i) = kin_ene(i) + pot_ene(i);
            Lz(i) = m_ear * (x(i)*vy(i) - y(i)*vx(i));
        //time arrows
            t(i) = t(i-1) + h;
//            if (r(i)<r(i-1)){
//                s = i;
//                jj = 1;
//                cout << "The Earth turns back after " << t(i)<< "years" <<endl;
//                emin = emid;
//                emid = (emax + emin)/2;
//                cout << "Inside if s=" << s<<endl;
//                cout <<fixed << setprecision(8)<<"Inside if emid=" << emid <<endl;
//                cout <<fixed << setprecision(8)<< "total energy= " << tot_ene(i) <<endl;
//                break;
//            }
        }
//        if (jj==0){
//            emax = emid;
//            emid = (emax + emin)/2;
//            cout <<fixed << setprecision(8)<<"Outside if emid=" <<  emid <<endl;
//            cout <<fixed << setprecision(8)<< "total energy= " << tot_ene(num_steps-1) <<endl;
//        }
//    }
    clock_t ends = clock();
    double elapsed_time = (double)(ends - start)/ CLOCKS_PER_SEC * 1000.0;
    cout << "Verlt"<<elapsed_time <<endl;
    //Write into file "Verlet_XX"
    Outputfile("Verlet", num_steps, elapsed_time, t_end, t, x, y, r, vx, vy, kin_ene, pot_ene, tot_ene, Lz);
}

// This function calculates time-evolution by Euler forward method.
void Euler(int num_steps, double t_end, double pi){
    double h = t_end /num_steps;                       //time steps
    double fac = 4.0 * h * pi * pi;                    //prefactor for acceleration
    double m_ear = 0.000003;                           //earth mass (in unit of M_sun(Sun mass) = 1)
    double G = 4*pi*pi;                                //Gravatational constant(in unit of AU^3*M_sun^-1*yr^-2 )
    double kin_fac = 0.5*m_ear, pot_fac = -G*m_ear;    //prefactors for kinetic and potential energy
    num_steps += 1;
    //Define and initialize coordinates x, y and their speed vx and vy
    vec x(num_steps), y(num_steps), vx(num_steps), vy(num_steps), t(num_steps), r(num_steps);
    x(0) = 1.0; y(0) = 0.0; vx(0) = 0.0;  vy(0)= 2.0*pi; t(0) = 0.0;
    r(0) = 1.0;
    //Define and initialize kinetic, potential and total energy (in unit of M_sun*AU^2*yr^-2),
    //as well as angular momentum Lz (in unit of M_sun*AU^2*yr^-1)
    vec kin_ene(num_steps), pot_ene(num_steps), tot_ene(num_steps), Lz(num_steps);
    kin_ene(0) = kin_fac*vy(0)*vy(0); pot_ene(0) = pot_fac/r(0); tot_ene(0) =  kin_ene(0)+pot_ene(0);
    Lz(0) = m_ear*vy(0);
    
    clock_t start = clock();
    
    for (int i = 1; i < num_steps; i++){
        //Calculate time-evolution of coordinates and velocities
        x(i) = x(i-1) + h * vx(i-1);
        y(i) = y(i-1) + h * vy(i-1);
        vx(i) = vx(i-1) - fac * x(i-1)/(r(i-1)*r(i-1)*r(i-1));
        vy(i) = vy(i-1) - fac * y(i-1)/(r(i-1)*r(i-1)*r(i-1));
        r(i) = sqrt(x(i) * x(i) + y(i) * y(i));
        //Calculate time-evolution of energies and angular momentum
        kin_ene(i) = kin_fac*(vx(i)*vx(i) + vy(i)*vy(i));
        pot_ene(i) = pot_fac/r(i);
        tot_ene(i) = kin_ene(i) + pot_ene(i);
        Lz(i) = m_ear * (x(i)*vy(i) - y(i)*vx(i));
        //time arrows
        t(i) = t(i-1) + h;
    }
    
    clock_t ends = clock();
    double elapsed_time = (double)(ends - start)/ CLOCKS_PER_SEC * 1000.0;
    cout << "Euler"<<elapsed_time <<endl;
    //Write into file "Euler_XX"
    Outputfile("Euler", num_steps, elapsed_time, t_end, t, x, y, r, vx, vy, kin_ene, pot_ene, tot_ene, Lz);
}

int main(int argc, char** argv)
{
    if (argc < 2){
        cout<<"More argument needed" <<endl;
        return 1;
    }
    
    double pi = acos(-1.0);
    double t_end = atof(argv[1]);      // time starts from 0, ends at t_end (in unit of year).
    int num_steps = 100000;              // number of time steps
    
    Verlet(num_steps, t_end, pi);
    Euler(num_steps, t_end, pi);
    
    return 0;
}
