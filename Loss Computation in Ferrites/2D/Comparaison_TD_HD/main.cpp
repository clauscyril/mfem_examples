#include <iostream>
#include <complex>
#include "mfem.hpp"

#include "../headers/customcurl.hpp" 
#include "../headers/Ferrite.hpp"
#include "../Harmonic_Domain/solver.hpp"
#include "../Time_Domain/solver_TD.hpp"

using namespace mfem;

int main() {
    const char *path = "../../../mesh/square.msh";
    Mesh *mesh = new Mesh(path, 1, 1); 

    // Defining the differents Ferrites and their parameters
    Ferrite N30("N30", 5.98e-2, 4.44e-1, 2.48e-6, 4300);
    Ferrite N87("N87", 4.24e-2, 1.48e-1, 2.68e-6, 2200);
    Ferrite T38("T38", 4.04e-2, 1.07e1, 8.06e-6, 10000);


    std::string name = "./data/HD_data.csv";                // Path to csv file for python plot
    std::ofstream data_file(name);                          // ofstream for writing in the file
    data_file << "f;p_eddy;flux_r;flux_i\n";                // Intialising the file with coluns names
    real_t P_eddy, P_mag;
    std::complex<real_t> flux = 0;
    real_t imax = 0;


    // Frequency range for the simulation
    real_t f = 100e3;
    real_t fc_end = 2e6;
    int N = 10; // N+1 points

    // Change of variable to obtain points linearly spaced on a log scale (optional)
    real_t u = log(f);
    real_t u_end = log(fc_end);
    real_t delta_u = (u_end - u)/ N; 

    for (int i = 0; i < N + 1; i++) {    // Looping the differents frequencies
        f = exp(u + i*delta_u);                            // Frequency for the simulation (from the change of variable)
        P_eddy = 0;
        P_mag = 0;
        imax = 0;
        GetPowerLoss(path, f, 1.8e6, P_eddy, P_mag, flux, imax);
        data_file << f << ";" << P_eddy << ";" << flux.real() << ";" << flux.imag() << std::endl;
        
        real_t I_rms = imax/sqrt(2);        // For now, the source term is a sine wave 
        real_t Ts = 1/f / 100;              // 100 echantillons par période
        int nb_period = 6;

        TD_sim(mesh, I_rms, f, Ts, (int)nb_period/(f*Ts), N30, 0);
    }
    return 0;
}