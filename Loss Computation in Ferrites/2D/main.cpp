#include <iostream>
#include <complex>
#include "mfem.hpp"

#include "solver.hpp"

#define M_PI 3.14159265358979323846

using namespace mfem;



// L'équation général est la suivante : 
//  rot(rho_eq rot(H)) + mu_eq jw H = 0
// Cependant, le problème présente une axisymétrie  => H = H_theta(r,z) . e_theta
// On peut donc se contenter d'un problème 2D 
// Malheureusment, lorsque H = H_theta(r,z) . e_theta, on n'a pas l'égalité suivante : 
// rot(rot H) . e_theta = - div(grad(H_theta)) mais on a rot(rot H) . e_theta = - div(grad(H_theta)) + H_theta/r²

// L'équation a résoudre est donc : 
//  - div(rho_eq grad(H_theta)) - rho_eq/r² * H_theta + jw mu_eq H_theta = 0


int main(){
    const char *path = "../../mesh/square.msh";


    for (int j = 0; j < 5; j++){
        std::string name = "./data";
        name += std::to_string(j) + ".csv";
        std::cout << name << std::endl;
        std::ofstream data_file(name);
        data_file << "fc;Ploss;flux\n";

        real_t fc_0 = 50e3;
        real_t fc = fc_0;
        real_t fc_end = 2e6;
        int N = 49;

        real_t fc_mu =  1e6;
        real_t fc_mu_end = 10e6;
        real_t delta_fc_mu = (fc_mu_end-fc_mu)/4;
        fc_mu = fc_mu + j*delta_fc_mu;
        real_t delta_f = (fc_end - fc)/N;
        std::cout << fc_mu << std::endl;

        real_t u = log(fc_0);
        real_t u_end = log(fc_end);
        real_t delta_u = (u_end - u)/ N;

        for (int i = 0; i < N + 1; i++) {
            fc = exp(u + i*delta_u);
            real_t PLoss, flux;
            GetPowerLoss(path, fc, fc_mu, PLoss, flux);
            data_file << fc << ";" << PLoss << ";" << flux <<std::endl; 
        }
    }

    // GetPowerLoss(path, 1e6, 5e6);
    // GetPowerLoss(path, 50e3, 5e6);
    return 0;
}

