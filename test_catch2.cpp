#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file
#include "entropy.h"
#include "catch2/catch.hpp"

TEST_CASE( "Entropy leaving empty cells is tested", "[entropy]") {
    std::vector<double>  t(4,0.0);
    REQUIRE(std::fabs(dropsize(t, 2)<=1e-8));
    double ent=0.693147;
    t[1]=0.5;
    t[0]=0.5;
    REQUIRE(std::fabs(dropsize(t,2)-0.5)<=1e-6);
    REQUIRE(std::fabs(Entropy(ent, 0, 1, 2, 3, t, 4, 2)-ent)<=1e-6);
}
