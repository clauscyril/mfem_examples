#include <iostream>
#include <complex>
#include "mfem.hpp"

#include "solver.hpp"

#define M_PI 3.14159265358979323846

int main(){
    const char *path = "../../mesh/square.msh";
    
    for (int j = 0; j < 5; j++){
        // Path to csv files for saving results
        std::string name = "./datafolder/data";
        name += std::to_string(j) + ".csv";

        // ofstream to those files
        std::ofstream data_file(name);
        data_file << "fc;Ploss;flux\n";

        // Frequency range for the simulation
        real_t fc_0 = 50e3;
        real_t fc = fc_0;
        real_t fc_end = 2e6;
        int N = 24; // N+1 points

        // Frequency range of the cutoff frequency of the permeability model
        real_t fc_mu =  1e6;
        real_t fc_mu_end = 10e6;
        real_t delta_fc_mu = (fc_mu_end-fc_mu)/4;
        fc_mu = fc_mu + j*delta_fc_mu;

        // Change of variable to obtain points linearly spaced on a logarithmic scale
        real_t u = log(fc_0);
        real_t u_end = log(fc_end);
        real_t delta_u = (u_end - u)/ N;

        for (int i = 0; i < N + 1; i++) {
            fc = exp(u + i*delta_u);                            // Frequency for the simulation
            real_t PLoss, P_mag;  
            std::complex<real_t> flux;
            real_t imax = 0;                               // parameters computed
            GetPowerLoss(path, fc, fc_mu, PLoss, P_mag, flux, imax);         // computation of the parameters
            data_file << fc << ";" << PLoss << ";" << flux <<std::endl; 
        }
    }
}