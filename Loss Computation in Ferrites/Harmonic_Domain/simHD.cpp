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
    std::string path = "../../mesh/square.msh";             // Path to the mesh
    std::string material = "N30";
    std::string b_option = "20";
    std::string f_value = "500e3";
    int viz = 0;


    OptionsParser args(argc, argv);
    args.AddOption(&path, "-p", "--path", "Path of the mesh file");
    args.AddOption(&f_value, "-f", "--Frequency", "Frequency of the input flux / current");
    args.AddOption(&material, "-mat", "--material", "Material of the ferrite");
    args.AddOption(&b_option, "-b", "--b_peak", "Average of b peak (mT)");
    args.AddOption(&viz, "-v", "--visualization", "0 : visualization disabled, 1 : enabled using GLVis");
    args.ParseCheck();

    
    Mesh *mesh = new Mesh(path, 1, 1);
    mesh->UniformRefinement();

    
    Ferrite N30("N30", 5.98e-2, 4.44e-1, 2.48e-6, 5111, 2277.94e3);
    Ferrite N87("N87", 4.24e-2, 1.48e-1, 2.68e-6, 2200, 1.8e6);
    Ferrite T38("T38", 4.04e-2, 1.07e1, 8.06e-6, 10000, 380e3);

    Ferrite ferrite = N30;

    real_t b_peak = 10e-3;

    if(material==N87.name) {
        ferrite = N87;
    } else if(material==T38.name) {
        ferrite = T38;
    } else if (material!=N30.name){
        std::cout <<  "Invalid ferrite name, using N30" << std::endl;
    }


    // Switch for selecting the value of b_peak used for the simulation
    if (b_option == "10")
    {
        b_peak = (real_t)10e-3;
    } else if (b_option == "20") {
        b_peak = (real_t)20e-3;
    } else {
        std::cout << "Invalid option for b_peak. Using b_peak = 10 mT" << std::endl;
    }
    
    real_t f = std::stod(f_value);
    real_t PLoss_eddy, PLoss_mag;
    std::complex<real_t> flux(0,0);    
    real_t NI=0;                            // parameters to be computed (will be passed as reference)
    GetPowerLoss(mesh, f, b_peak, NI, ferrite, PLoss_eddy, PLoss_mag, flux, viz);   // Calling the functions that computes the power losses
    std::cout << f << ";" << PLoss_eddy << ";" << PLoss_mag << ";" << PLoss_eddy + PLoss_mag << ";" << flux.real() << ";" << flux.imag() << ";"<< NI << std::endl;

return 0;
}

