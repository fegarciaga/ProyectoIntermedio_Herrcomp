#include <iostream>
#include <cstdlib>
#include <string>
#include "entropy.h"

int main(int argc, char* argv[])
{
    int N_mol,  Nbin, seed1;
    double l, t;
    std::ifstream fin("parameters.txt");
    std::string tmp;
    fin >> N_mol >> tmp;
    fin >> l >> tmp;
    fin >> t >> tmp;
    fin >> seed1 >> tmp;
    fin.close();
    std::cout<< t <<"\t"<<N_mol << "\t"<<seed1 <<"\t" << Nbin<<"\t"<<l<<"\n";
    Nbin=100;
    simulation(t, N_mol, seed1, Nbin, l);
    simulationwithhole(t, N_mol, seed1, Nbin, l);
    return 0;
}
