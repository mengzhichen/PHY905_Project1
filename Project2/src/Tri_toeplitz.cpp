#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include "Outputfile.hpp"
#include "armadillo"

using namespace std;
using namespace arma;

void maxoffele(mat A, int *r, int *c, int n){
    max = 0;
    for (int i = 0; i < n; i++){
        for (int j = i+1; j < n; j++){
            max_temp = A(i,j) * A(i,j);
            if (max_temp > max){
                max = max_temp;
                r = i;
                c = j;
            }
        }
}
                            
void jocabi(mat A, mat B, int r, int c, int n){
    double tau, tan, cos2, sincos, sin2, cos, sin;
    tau = (A(c,c) - A(r,r))/2/A(r,c));               //cot(2theta)
    if (tau >= 0) tan = -tau + sqrt(1 + tau * tau);
    else tan = -tau - sqrt(1 + tau * tau);
    tan = -tau - sqrt(1 + tau * tau);
    cos2 = 1/(1 + tan * tan);
    sincos = cos2 * tan;
    sin2 = sincos * tan;
    cos = 1/(sqrt(1 + tan*tan));
    sin = tan * cos;
    B(c,c) = A(c,c)*cos2 - 2*A(c,r)*sincos + A(r,r)*sin2;
    B(r,r) = A(r,r)*cos2 + 2*A(c,r)*sincos + A(c,c)*sin2;
    B(r,c) = (A(c,c)-A(r,r))*sincos + A(c,r)*(cos2 - sin2);
    for (int i = 0; i < n; i++){
        if (i != r && i != c){
            B(i,c) = A(i,c)*cos - A(i,r)*sin;
            B(i,r) = A(i,r)*cos + A(i,c)*sin;
        }
    }
                            
 void Cyclic_jacobi()
                            
                            
int main(int argc, char** argv)
{
    int exponent = 0;
    if (argc < 2){
        cout<<"More argument needed";
        return 1;
    }
    //string filename = "LU_decomposition";
    
    for(int i = 1; i < argc; i++){
        double tol = 1;
        
        
        clock_t start = clock();
        
;
    
        clock_t ends = clock();
        double elapsed_time = (double)(ends - start)/ CLOCKS_PER_SEC * 1000.0;
        
        //Outputfile("LU_soulution", n, elapsed_time, u, x);
       // delete [] u; delete [] y; delete [] y_temp; delete [] x;
    }
    return 0;
}
