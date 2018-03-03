#include "catch.hpp"
#include "ArmaDiag.h"

TEST_CASE("Test for eigenvalues by comparing with analytic solution"){
    int n = 8;
    vec eigval(n), ana_eigval(n);
    // calculate eigenvalues using Armadilo
    eigval = Toeplitz(n);
    // extract analytical eigenvalues
    ana_eigval = Analytical(n);
    REQUIRE(eigval(0)==Approx(ana_eigval(0)).epsilon(0.00001));
    REQUIRE(eigval(1)==Approx(ana_eigval(1)).epsilon(0.00001));
    REQUIRE(eigval(2)==Approx(ana_eigval(2)).epsilon(0.00001));
    REQUIRE(eigval(3)==Approx(ana_eigval(3)).epsilon(0.00001));
}
