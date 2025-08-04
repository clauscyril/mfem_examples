#include <iostream>
#include <complex>
#include "mfem.hpp"

#include "../headers/customcurl.hpp"
#include "solver.hpp"

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif

using namespace mfem;


int main(int argc, char *argv[]){
    const char *path = "../../mesh/square.msh";             // Path to the mesh

    Mesh *mesh = new Mesh(path, 1, 1);
    mesh->UniformRefinement();

    std::string material = "N30";
    Ferrite N30("N30", 5.98e-2, 4.44e-1, 2.48e-6, 4300, 1.8e6);
    Ferrite N87("N87", 4.24e-2, 1.48e-1, 2.68e-6, 2200, 1.8e6);
    Ferrite T38("T38", 4.04e-2, 1.07e1, 8.06e-6, 10000, 380e3);

    Ferrite ferrite = N30;

    std::string b_option = "10";
    real_t b_peak = 10e-3;

    OptionsParser args(argc, argv);
    args.AddOption(&material, "-m", "--material", "Material to use.");
    args.AddOption(&b_option, "-b", "--b_option", "10 or 20 mT");
    args.ParseCheck();


    if (material == N30.name){
        std::cout << "Using N30 ferrite." << std::endl;
    } else if(material==N87.name) {
        ferrite = N87;
        std::cout << "Using N87 ferrite." << std::endl;
    } else if(material==T38.name) {
        ferrite = T38;
        std::cout << "Using T38 ferrite." << std::endl;
    } else {
        std::cout <<  "Invalid ferrite name, using N30" << std::endl;
    }


    // Switch for selecting the value of b_peak used for the simulation
    if (b_option == "10")
    {
        b_peak = (real_t)10e-3;
        std::cout << "b_peak = 10 mT" << std::endl;
    } else if (b_option == "20") {
        b_peak = (real_t)20e-3;
    } else {
        std::cout << "Invalid option for b_peak. Using b_peak = 10 mT" << std::endl;
    }
    


    // ********************************
    // // Path to csv files for saving results
    std::string name = "../data/MFEM/unipg_mfem";   // Path to csv file for python plot
    name += ferrite.name + "_" +b_option + ".csv";
    std::ofstream data_file(name);                          // ofstream for writing in the file
    data_file << "fc;P_eddy;P_mag;P_tot;flux_r;flux_i;NI\n";                    // Intialising the file with coluns names

    
    // Frequency range for the simulation
    real_t f_0 = 100e3;
    real_t f = f_0;
    real_t f_end = 2e6;
    int N = 14; // N+1 points

    // Change of variable to obtain points linearly spaced on a log scale (optional)
    real_t u = log(f_0);
    real_t u_end = log(f_end);
    real_t delta_u = (u_end - u)/ N; 

    // real_t fc_mu = 2.5e6;  // cutoff frequency of µeq = µ0 µr / (1 + jw/wc) with wc = 2 pi fc

    for (int i = 0; i < N + 1; i++) {    // Looping the differents frequencies
        f = exp(u + i*delta_u);          // Frequency for the simulation (from the change of variable)
        real_t PLoss_eddy, PLoss_mag;
        std::complex<real_t> flux(0,0);    
        real_t NI=0;                            // parameters to be computed (will be passed as reference)
        GetPowerLoss(mesh, f, b_peak, NI, ferrite, PLoss_eddy, PLoss_mag, flux, false);   // Calling the functions that computes the power losses
        data_file << f << ";" << PLoss_eddy << ";" << PLoss_mag << ";" << PLoss_eddy + PLoss_mag << ";" << flux.real() << ";" << flux.imag() << ";"<< NI << std::endl; // Writing the results in the csv file
    }

return 0;
}

