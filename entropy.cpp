#include "entropy.h"

void simulation (double T_MAX, int N, int seed, int NBINS, double Nsize, std::string filename)
{

    std::ofstream fout(filename, std::ofstream::out);
    fout.precision(15);
    fout.setf(std::ios::scientific);
    
    //Creates random distributon
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    
    //Defines a vector with position
    std::vector<double> z(N*2,0.0);
    //Defines the lattice vector
    std::vector<double> lattice(NBINS*NBINS,0.0);
    double DX=2*Nsize/NBINS;
    double entr=0;

    //Fill vector with random position
    for (int ii=0; ii < N*2; ++ii) z[ii] = dis(gen);

    //Calculates initial probability density
    prob(lattice, z, N, NBINS, Nsize);
     
    double r=dropsize(z, N);
    
    //Calculates initial entropy
    for (int ii=0; ii<NBINS;++ii)
    {
        for (int jj=0; jj<NBINS; ++jj)
        {
            if(std::fabs(lattice[ii*NBINS+jj])>1e-10)
            {
                entr-=lattice[ii*NBINS+jj]*std::log(lattice[ii*NBINS+jj]);
            }
        }
    }

    fout << 0 << "\t"<<entr<<"\t"<< r<< "\n";
    //evolution
    for (int ii=1; ii<T_MAX; ++ii)
    {
        double aux1=0, aux2=0;
        int bin1=0, bin2=0, newbin1=0, newbin2=0;
        int particle=0;
        // choose a random particle
        particle=std::fabs(N*dis(gen));
        bin1=int((z[2*particle]+Nsize)/DX);
        bin2=int((z[2*particle+1]+Nsize)/DX);
        // check that the particle stays with in the boundary
        aux1=z[particle*2]+dis(gen);
        if (std::fabs(aux1)<Nsize)
        {
            aux2=z[particle*2+1]+dis(gen);
            if (std::fabs(aux2)<Nsize)
            {
                newbin1=int((aux1+Nsize)/DX);
                newbin2=int((aux2+Nsize)/DX);
                r*=r;
                r-=(z[particle*2]*z[particle*2]+z[particle*2+1]*z[particle*2+1])/N;
                r+=(aux1*aux1+aux2*aux2)/N;
                r=std::sqrt(r);
                z[particle*2]=aux1;
                z[particle*2+1]=aux2;
                //Calculates the entropy
                entr=Entropy(entr, bin1, bin2, newbin1, newbin2, lattice, NBINS, N);
            }
        }
        fout<<ii<<"\t"<<entr<<"\t"<<r<<"\n";
    }
    fout.close();
    
    return;
}

double Entropy (double initialvalue, int NBIN1, int NBIN2, int NNEWBIN1, int NNEWBIN2, std::vector<double> &distribution, int Ntotal, int size)
{
    if(std::fabs(distribution[NNEWBIN1*Ntotal+NNEWBIN2])>1e-10)
    {
        initialvalue+=distribution[NBIN1*Ntotal+NBIN2]*std::log(distribution[NBIN1*Ntotal+NBIN2])+distribution[NNEWBIN1*Ntotal+NNEWBIN2]*std::log(distribution[NNEWBIN1*Ntotal+NNEWBIN2]);
        distribution[NBIN1*Ntotal+NBIN2]-=1.0/size;
        distribution[NNEWBIN1*Ntotal+NNEWBIN2]+=1.0/size;
        if(std::fabs(distribution[NBIN1*Ntotal+NBIN2])>1e-10)
        {
            initialvalue-=distribution[NBIN1*Ntotal+NBIN2]*std::log(distribution[NBIN1*Ntotal+NBIN2])+distribution[NNEWBIN1*Ntotal+NNEWBIN2]*std::log(distribution[NNEWBIN1*Ntotal+NNEWBIN2]);
        }
        else
        {
            initialvalue-=distribution[NNEWBIN1*Ntotal+NNEWBIN2]*std::log(distribution[NNEWBIN1*Ntotal+NNEWBIN2]);
        }
    }
    else
    {
        initialvalue+=distribution[NBIN1*Ntotal+NBIN2]*std::log(distribution[NBIN1*Ntotal+NBIN2]);
        distribution[NBIN1*Ntotal+NBIN2]-=1.0/size;
        distribution[NNEWBIN1*Ntotal+NNEWBIN2]+=1.0/size;
        if(std::fabs(distribution[NBIN1*Ntotal+NBIN2])>1e-10)
        {
            initialvalue-=distribution[NBIN1*Ntotal+NBIN2]*std::log(distribution[NBIN1*Ntotal+NBIN2])+distribution[NNEWBIN1*Ntotal+NNEWBIN2]*std::log(distribution[NNEWBIN1*Ntotal+NNEWBIN2]);
        }
        else
        {
            initialvalue-=distribution[NNEWBIN1*Ntotal+NNEWBIN2]*std::log(distribution[NNEWBIN1*Ntotal+NNEWBIN2]);
        }
    }
    return initialvalue;
}

double dropsize (std::vector<double> & position, int size)
{
    double distance=0;
    for (int ii=0; ii<2*size; ++ii)
    {
        distance += position[ii]*position[ii];
    }
    return std::sqrt(distance/size);
}

void prob (std::vector<double> &distribution, std::vector<double> &position, int psize, int dsize, double range)
{
    double d=2*range/dsize;
    for (int ii=0; ii<psize; ++ii)
    {
        int bin1=int((position[2*ii]+range)/d);
        int bin2=int((position[2*ii+1]+range)/d);
        distribution[bin1*dsize+bin2]+=1.0/psize;
    }

    return;
}

void simulationwithhole (double T_MAX, int N, int seed, int NBINS, double Nsize, std::string filename)
{

    std::ofstream fout(filename, std::ofstream::out);
    fout.precision(15);
    fout.setf(std::ios::scientific);
    
    //Creates random distributon
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    
    //Defines a vector with position
    std::vector<double>z(N*2, 0.0);
    //Defines the lattice vector
    std::vector<double> lattice(NBINS*NBINS+1,0.0);
    double DX=2*Nsize/NBINS;

    //Fill vector with random position
    for (int ii=0; ii < N*2; ++ii) z[ii] = dis(gen);

    //Calculates initial probability density
    prob(lattice, z, N, NBINS, Nsize);
    int n=N;
    fout<< 0<< "\t"<<n<<"\n";

    //evolution
    for (int ii=1; ii<T_MAX; ++ii)
    {
        double aux1=0, aux2=0;
        int bin1=0, bin2=0, newbin1=0, newbin2=0;
        int particle=0;
        // choose a random particle
        particle=std::fabs(N*dis(gen));
        bin1=int((z[2*particle]+Nsize)/DX);
        bin2=int((z[2*particle+1]+Nsize)/DX);
        //Applies algorithm only if the particle is inside the box
        if(std::fabs(z[2*particle])<Nsize)
        {
            // check that the particle stays with in the boundary
            aux1=z[particle*2]+dis(gen);
            aux2=z[particle*2+1]+dis(gen);
            if (std::fabs(aux1)<Nsize)
            {
                if (std::fabs(aux2)<Nsize)
                {
                    newbin1=int((aux1+Nsize)/DX);
                    newbin2=int((aux2+Nsize)/DX);
                    z[particle*2]=aux1;
                    z[particle*2+1]=aux2;
                    lattice[bin1*NBINS+bin2]-=1.0/N;
                    lattice[newbin1*NBINS+newbin2]+=1.0/N;

                }
            }
            if(aux1>Nsize)
            {
                if(std::fabs(aux2)<Nsize/10)
                {
                    z[particle*2]=aux1;
                    z[particle*2+1]=aux2;
                    lattice[bin1*NBINS+bin2]-=1.0/N;
                    lattice[NBINS*NBINS]+=1.0/N;
                    n-=1;
                }
            }
        }
        fout << ii<<"\t"<<n<<"\n";
    }
    fout.close();
    return;
}

