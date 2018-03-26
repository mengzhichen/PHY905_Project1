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
void Outputfile(string, int, double, vec, vec, vec, vec, vec, vec, vec, vec, vec);
void Verlet( int, double, double);
void Euler( int, double, double);

// This function write results into file.
void Outputfile(string filename, int num_steps, double t_end, vec t, vec x, vec y, vec vx, vec vy, vec kin, vec pot, vec tot, vec Lz){
    string fileout = filename;
    string size = to_string(t_end);
    fileout.append("_"+size+".txt");
    ofstream file(fileout);
    file << "time(yr) " << setw(5) << "x" << setw(10) << "y" << setw(10) << "vx"<< setw(10) << "vy"  << setw(13) << "kin_ene" << setw(10) << "pot_ene"<< setw(10) << "tot_ene" << setw(8) << "Lz"<< endl;
    cout << setiosflags(ios::fixed);
    for (int i = 0; i < num_steps; i++){
        file << fixed << setprecision(6) << t(i) << setw(10) << x(i) << setw(10) << y(i) << setw(10) << vx(i) << setw(10) << vy(i) << setw(10) << kin(i) << setw(10) << pot(i) << setw(10) << tot(i) << setw(10) << Lz(i) << endl;
    }
    file.close();
}

// This function calculates time-evolution by velovity Verlet method.
void Verlet(int num_steps, double t_end, double pi){
    double h = t_end /num_steps;                       //time steps
    double fac = 4.0 * h * pi * pi;                    //prefactor for acceleration
    double m_ear = 0.000003;                           //Earth mass (in unit of M_sun(Sun mass) = 1)
    double mj = 0.00095;                               //Jupiter mass (in unit of M_sun(Sun mass) = 1)
    double rj = 5.20;                                  //Sun Jupiter distance (in unit of AU)
    double G = 4*pi*pi;                                //Gravatational constant(in unit of AU^3*M_sun^-1*yr^-2 )
    double kin_fac = 0.5*m_ear, pot_fac = -G*m_ear;    //prefactors for kinetic and potential energy
    num_steps += 1;
    //Define and initialize coordinates x, y and their speed vx and vy for earth
    vec x(num_steps), y(num_steps), vx(num_steps), vy(num_steps), t(num_steps);
    x(0) = 1.0; y(0) = 0.0; vx(0) = 0.0;  vy(0)= 2.0*pi; t(0) = 0.0;
    double r = 1.0;
    //Define and initialize coordinates x, y and their speed vx and vy for jupiter
    vec xj(num_steps), yj(num_steps), vxj(num_steps), vyj(num_steps);
    xj(0) = 5.2; yj(0) = 0.0; vxj(0) = 0.0;  vyj(0)= 2.0*pi/(sqrt(rj));
    double xej = 4.2, yej = 0.0, rej = 4.2;
    //Define and initialize kinetic, potential and total energy (in unit of M_sun*AU^2*yr^-2),
    //as well as angular momentum Lz (in unit of M_sun*AU^2*yr^-1)
    vec kin_ene(num_steps), pot_ene(num_steps), tot_ene(num_steps), Lz(num_steps);
    kin_ene(0) = kin_fac*vy(0)*vy(0); pot_ene(0) = pot_fac/r; tot_ene(0) =  kin_ene(0)+pot_ene(0);
    Lz(0) = m_ear*vy(0);
    //acceleration times h at step i and i+1 for earth;
    double axh0 = 0.0, ayh0 = 0.0, axh1 = 0.0, ayh1 = 0.0;
    //acceleration times h at step i and i+1 for jupiter;
    double axhj0 = 0.0, ayhj0 = 0.0, axhj1 = 0.0, ayhj1 = 0.0;
    
    for (int i = 1; i < num_steps; i++){
        //Calculate time-evolution of coordinates and velocities
        //Calculate distances
        r = sqrt(x(i-1) * x(i-1) + y(i-1) * y(i-1));
        rj = sqrt(xj(i-1) * xj(i-1) + yj(i-1) * yj(i-1));
        xej = xj(i-1) - x(i-1);
        yej = yj(i-1) - y(i-1);
        rej = sqrt(xej*xej + yej*yej);
        //Calculate coordinates and acceleration factors
        axh0 = fac * ( x(i-1)/(r*r*r) - mj*xej/(rej*rej*rej) );
        ayh0 = fac * ( y(i-1)/(r*r*r) - mj*yej/(rej*rej*rej) );
        axhj0 = fac * ( xj(i-1)/(rj*rj*rj) + m_ear*xej/(rej*rej*rej) );
        ayhj0 = fac * ( yj(i-1)/(rj*rj*rj) + m_ear*yej/(rej*rej*rej) );
        
        x(i) = x(i-1) + h * vx(i-1) - h * axh0/2;
        y(i) = y(i-1) + h * vy(i-1) - h * ayh0/2;
        xj(i) = xj(i-1) + h * vxj(i-1) - h * axhj0/2;
        yj(i) = yj(i-1) + h * vyj(i-1) - h * ayhj0/2;
        
        r = sqrt(x(i) * x(i) + y(i) * y(i));
        rj = sqrt(xj(i) * xj(i) + yj(i) * yj(i));
        xej = xj(i) - x(i);
        yej = yj(i) - y(i);
        rej = sqrt(xej*xej + yej*yej);
        
        axh1 = fac * ( x(i)/(r*r*r) - mj*xej/(rej*rej*rej) );
        ayh1 = fac * ( y(i)/(r*r*r) - mj*yej/(rej*rej*rej) );
        axhj1 = fac * ( xj(i)/(rj*rj*rj) + m_ear*xej/(rej*rej*rej) );
        ayhj1 = fac * ( yj(i)/(rj*rj*rj) + m_ear*yej/(rej*rej*rej) );
        
        vx(i) = vx(i-1) - 0.5*(axh0+axh1);
        vy(i) = vy(i-1) - 0.5*(ayh0+ayh1);
        vxj(i) = vxj(i-1) - 0.5*(axhj0+axhj1);
        vyj(i) = vyj(i-1) - 0.5*(ayhj0+ayhj1);
        //Calculate time-evolution of energies and angular momentum
        kin_ene(i) = kin_fac*(vx(i)*vx(i) + vy(i)*vy(i));
        pot_ene(i) = pot_fac/r;
        tot_ene(i) = kin_ene(i) + pot_ene(i);
        Lz(i) = m_ear * (x(i)*vy(i) - y(i)*vx(i));
        //time arrows
        t(i) = t(i-1) + h;
    }
    //Write into file "Verlet_XX"
    Outputfile("SEJ_Verlet", num_steps, t_end, t, x, y, vx, vy, kin_ene, pot_ene, tot_ene, Lz);
}

int main(int argc, char** argv)
{
    if (argc < 2){
        cout<<"More argument needed" <<endl;
        return 1;
    }
    
    double pi = acos(-1.0);
    double t_end = atof(argv[1]);      // time starts from 0, ends at t_end (in unit of year).
    int num_steps = 5000;              // number of time steps
   
    //Euler(num_steps, t_end, pi);
    Verlet(num_steps, t_end, pi);
    return 0;
}
