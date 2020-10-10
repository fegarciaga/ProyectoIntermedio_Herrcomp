#include "entropy.h"

void simulation (double T_MAX, int N, int seed, int NBINS)
{
    //Creates random distributon
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    
    //Defines a vector with position
    double *z = new double [N*2] {0.0};

    //Defines the lattice vector
    double *lattice = new double [NBINS*NBINS] {0.0};
    double DX=20.0/NBINS;
    double entr=0;

    //Fill vector with random position
    for (int ii=0; ii < N*2; ++ii) z[ii] = dis(gen);

    //Calculates initial probability density
    for (int ii=0; ii<N; ++ii)
    {
        int bin1=int((z[2*ii]+10.0)/DX);
        int bin2=int((z[2*ii+1]+10.0)/DX);
        lattice[bin1*NBINS+bin2]+=1.0/N;
    }

    //Calculates initial density
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

    std::cout << 0 << "\t"<<entr<<"\n";
    //evolution
    for (int ii=1; ii<T_MAX; ++ii)
    {
        double aux1=0, aux2=0;
        int bin1=0, bin2=0, newbin1=0, newbin2=0;
        int particle=0;
        // choose a random particle
        particle=std::fabs(N*dis(gen));
        bin1=int((z[2*particle]+10.0)/DX);
        bin2=int((z[2*particle+1]+10.0)/DX);
        // check that the particle stays with in the boundary
        aux1=z[particle*2]+dis(gen);
        if (std::fabs(aux1)<10)
        {
            aux2=z[particle*2+1]+dis(gen);
            if (std::fabs(aux2)<10)
            {
                newbin1=int((aux1+10.0)/DX);
                newbin2=int((aux2+10.0)/DX);
                z[particle*2]=aux1;
                z[particle*2+1]=aux2;
                if(std::fabs(lattice[newbin1*NBINS+newbin2])>1e-10)
                {
                    entr+=lattice[bin1*NBINS+bin2]*std::log(lattice[bin1*NBINS+bin2])+lattice[newbin1*NBINS+newbin2]*std::log(lattice[newbin1*NBINS+newbin2]);
                    lattice[bin1*NBINS+bin2]-=1.0/N;
                    lattice[newbin1*NBINS+newbin2]+=1.0/N;
                    if(std::fabs(lattice[bin1*NBINS+bin2])>1e-10)
                    {
                        entr-=lattice[bin1*NBINS+bin2]*std::log(lattice[bin1*NBINS+bin2])+lattice[newbin1*NBINS+newbin2]*std::log(lattice[newbin1*NBINS+newbin2]);
                    }
                    else
                    {
                        entr-=lattice[newbin1*NBINS+newbin2]*std::log(lattice[newbin1*NBINS+newbin2]);
                    }
                }
                else
                {
                    entr+=lattice[bin1*NBINS+bin2]*std::log(lattice[bin1*NBINS+bin2]);
                    lattice[bin1*NBINS+bin2]-=1.0/N;
                    lattice[newbin1*NBINS+newbin2]+=1.0/N;
                    if(std::fabs(lattice[bin1*NBINS+bin2])>1e-10)
                    {
                        entr-=lattice[bin1*NBINS+bin2]*std::log(lattice[bin1*NBINS+bin2])+lattice[newbin1*NBINS+newbin2]*std::log(lattice[newbin1*NBINS+newbin2]);
                    }
                    else
                    {
                        entr-=lattice[newbin1*NBINS+newbin2]*std::log(lattice[newbin1*NBINS+newbin2]);
                    }
                }
            }
        }
        //std::cout<<ii<<"\t"<<entr<<"\n";
        if(ii==7e5)
        {
            for (int jj=0; jj<N; ++jj)
            {
                std::cout<<z[2*jj]<<"\t"<<z[2*jj+1]<<"\n";
            }
        }
    }
    delete [] z;
    delete [] lattice;
    
    return;
}
