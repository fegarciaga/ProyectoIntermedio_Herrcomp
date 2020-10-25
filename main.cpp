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
    Nbin=100;
    std::string sim="datossimulation.txt", hol="datoshole.txt";
    simulation(t, N_mol, seed1, Nbin, l, sim);
    simulationwithhole(t, N_mol, seed1, Nbin, l, hol);
    /* loop for finding the equilibrium time
    for (int ii=5; ii<=25; ii+=5)
    {
        std::string fname = "datost_"+std::to_string(ii)+".txt";
        simulation(t, N_mol, seed1, Nbin, ii, fname);
        }*/
    return 0;
}
