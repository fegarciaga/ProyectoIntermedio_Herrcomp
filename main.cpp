#include <iostream>
#include <cstdlib>
#include "entropy.h"

int main(int argc, char* argv[])
{
    double t=std::atof(argv[1]);
    int N_mol=std::atoi(argv[2]);
    int seed1=std::atoi(argv[3]);
    int Nbin=std::atoi(argv[4]);
    double l=std::atof(argv[5]);
    simulation(t, N_mol, seed1, Nbin, l);
    return 0;
}
