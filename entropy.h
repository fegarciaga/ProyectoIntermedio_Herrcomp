#include <iostream>
#include <random>
#include <cmath>

void simulation (double T_MAX, int N, int seed, int NBINS, double Nsize);
double entropy (double initialvalue, int NBIN1, int NBIN2, int NNEWBIN1, int NNEWBIN2, double *distribution, int Ntotal, int size);
double dropsize (double * position, int Nsize);
void prob (double *distribution, double *position, int psize, int dsize, double range);
void simulationwithhole (double T_MAX, int N, int seed, int NBINS, double Nsize);
int Number (double * distribution, int dsize, int psize);
