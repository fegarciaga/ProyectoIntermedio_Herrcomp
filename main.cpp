#include <iostream>
#include <cstdlib>
#include <string>
#include <chrono>
#include "entropy.h"

int main(int argc, char* argv[])
{
    int N_mol,  Nbin, seed1;
    double l, t, time=0;
    std::ifstream fin("parameters.txt");
    std::string tmp;
    fin >> N_mol >> tmp;
    fin >> l >> tmp;
    fin >> t >> tmp;
    fin >> seed1 >> tmp;
    fin.close();
    Nbin=100;
    l=std::atof(argv[1]);
    std::string sim="datossimulation.txt", hol="datoshole.txt";
    int REPS = std::atoi(argv[2]);
    for(int ii=0; ii<REPS; ++ii)
    {
        auto start = std::chrono::steady_clock::now();
        simulation(t, N_mol, seed1, Nbin, l, sim);
        simulationwithhole(t, N_mol, seed1, Nbin, l, hol);
        auto end =std::chrono::steady_clock::now();
        time+=std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()*1.0e-3;
    }
    /* loop for finding the equilibrium time
    for (int ii=5; ii<=25; ii+=5)
    {
        std::string fname = "datost_"+std::to_string(ii)+".txt";
        simulation(t, N_mol, seed1, Nbin, ii, fname);
        }*/
    std::cout  << time/REPS<<"\n";
    return 0;
}
