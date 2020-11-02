#ifndef __ENTROPY_H_
#define __ENTROPY_H_

#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

void simulation (double T_MAX, int N, int seed, int NBINS, double Nsize, std::string filename);
double Entropy (double initialvalue, int NBIN1, int NBIN2, int NNEWBIN1, int NNEWBIN2, std::vector<double> &distribution, int Ntotal, int size);
double dropsize (std::vector<double> & position, int Nsize);
void prob (std::vector<double> &distribution, std::vector<double> &position, int psize, int dsize, double range);
void simulationwithhole (double T_MAX, int N, int seed, int NBINS, double Nsize, std::string filename);

#endif
