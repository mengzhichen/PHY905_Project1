#include "catch.hpp"
#include "ClassicJacobi.h"


TEST_CASE("Testing max off-diagonal element a(i,j)"){
    int n = 5;
    double offmax = 0.0, norm = 0.0;
    int r = 0, c = 0;
    mat A = zeros<mat>(n-1,n-1);
    //initialize matrices and vector
    Initialization(A, norm, n);
    //find maximum off-diagonal matrix element
    maxoffele(A, r, c, offmax, n);
    
    REQUIRE(r==0);
    REQUIRE(c==1);
}

TEST_CASE("Test the conservation of Frobenius norm"){
    int n = 5, iter = 0, maxiter = 100;
    double offmax = 10.0, norm = 0.0, norm2 = 0.0;
    int r = 0, c = 0;
    double tol = 0.001;
    mat A = zeros<mat>(n-1,n-1);
    mat eigvec = eye<mat>(n-1,n-1);
    //initialize matrix; set tolerence; 
    Initialization(A, norm, n);
    //do Jacobi iteration until convergence
    while (offmax > tol && iter < maxiter){
            maxoffele(A, r, c, offmax, n);
            jacobi(A, eigvec, r, c, n);
            iter += 1;
    }
    Frobeniusnorm(A, norm2, n);
    double diff = norm - norm2;
    REQUIRE(diff == Approx(0.0).epsilon(0.001));
}

TEST_CASE("Test for eigenvalues; correct eigenvectors and orthonormality"){
    int n = 4, iter = 0, maxiter = 100;
    double offmax = 10.0, norm = 0.0, norm2 = 0.0;
    int r = 0, c = 0;
    double tol = 0.001;
    mat A = zeros<mat>(n-1,n-1);
    mat eigvec = eye<mat>(n-1,n-1);
    vec eigval(n-1);
    //initialize matrices and vector
    Initialization(A, norm, n);
    //do Jacobi iteration until convergence
    while (offmax > tol && iter < maxiter){
            maxoffele(A, r, c, offmax, n);
            jacobi(A, eigvec, r, c, n);
            iter += 1;
    }
    ///// place eigenvalues and eigenvecotrs in sequence
    for (int r = 0; r < n - 1; r++){
        for (int c = r; c < n - 1; c++){
            if (A(r,r) > A(c,c)){
                double A_temp = A(c,c);
                A(c,c) = A(r,r);
                A(r,r) = A_temp;
                eigvec.swap_cols(r,c);
             }
        }
    }
    // extract eigenvalues
    for (int i = 0; i < n-1; i++){
            eigval(i) = A(i,i);
    }
    // test for eigenvalues
    REQUIRE(eigval[0]==Approx(0.5858));
    REQUIRE(eigval[1]==Approx(2.0000));
    REQUIRE(eigval[2]==Approx(3.4142));

    // test for correctness of eigenvalues
    REQUIRE(eigvec(0,0)==Approx(0.50).epsilon(0.01));
    REQUIRE(eigvec(0,1)==Approx(-0.70711).epsilon(0.01));
    REQUIRE(eigvec(0,2)==Approx(0.50).epsilon(0.01));
    REQUIRE(eigvec(1,0)==Approx(0.70711).epsilon(0.01));
    REQUIRE(eigvec(1,1)==Approx(0.0).epsilon(0.01));
    REQUIRE(eigvec(1,2)==Approx(-0.70711).epsilon(0.01));
    REQUIRE(eigvec(2,0)==Approx(0.50).epsilon(0.01));
    REQUIRE(eigvec(2,1)==Approx(0.70711).epsilon(0.01));
    REQUIRE(eigvec(2,2)==Approx(0.50).epsilon(0.01));

    // test for orthonormality of eigenvalues
    //dot0=v0*v1=0
    double dot0=eigvec(0,0)*eigvec(1,0)+eigvec(0,1)*eigvec(1,1)
        +eigvec(0,2)*eigvec(1,2);
    //dot1=v0*v2=0
    double dot1=eigvec(0,0)*eigvec(2,0)+eigvec(0,1)*eigvec(2,1)
        +eigvec(0,2)*eigvec(2,2);
    //dot2=v1*v2=0
    double dot2=eigvec(1,0)*eigvec(2,0)+eigvec(1,1)*eigvec(2,1)
        +eigvec(1,2)*eigvec(2,2);
    //vecnorm0=v0*v0=1
    double vecnorm0=eigvec(0,0)*eigvec(0,0)+eigvec(0,1)*eigvec(0,1)
        +eigvec(0,2)*eigvec(0,2);
    //vecnorm1=v1*v1=1
    double vecnorm1=eigvec(1,0)*eigvec(1,0)+eigvec(1,1)*eigvec(1,1)
        +eigvec(1,2)*eigvec(1,2);
    //vecnorm2=v2*v2=1
    double vecnorm2=eigvec(2,0)*eigvec(2,0)+eigvec(2,1)*eigvec(2,1)
        +eigvec(2,2)*eigvec(2,2);
    REQUIRE(dot0==Approx(0.000).epsilon(0.01));
    REQUIRE(dot1==Approx(0.000).epsilon(0.01));
    REQUIRE(dot2==Approx(0.000).epsilon(0.01));
    REQUIRE(vecnorm0==Approx(1.000).epsilon(0.01));
    REQUIRE(vecnorm1==Approx(1.000).epsilon(0.01));
    REQUIRE(vecnorm2==Approx(1.000).epsilon(0.01));
}
