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
// rot(rot H) . e_theta = - div(grad(H_theta)) mais on a rot(rot H) . e_theta = - div(grad(H_theta)) - H_theta/r²

// L'équation a résoudre est donc : 
//  - div(rho_eq grad(H_theta)) - rho_eq/r² * H_theta + jw mu_eq H_theta = 0


int main(){
    const char *path = "../../mesh/square.msh";
    for (int j = 0; j < 5; j++){
        std::string name = "./data";
        name += std::to_string(j) + ".csv";
        std::cout << name << std::endl;
        std::ofstream data_file(name);
        data_file << "fc;Ploss\n";
        real_t fc = 100e3;
        real_t fc_end = 2e6;
        int N = 9;
        real_t delta_fc_mu = (10e6-1e6)/4;
        real_t fc_mu = 1e6 + j*delta_fc_mu;
        real_t delta_f = (fc_end - fc)/N;
        std::cout << fc_mu << std::endl;
        for (int i = 0; i < N + 1; i++) {
            real_t PLoss = GetPowerLoss(path, fc, fc_mu);
            data_file << fc << ";" << PLoss << std::endl; 
            fc += delta_f;
        }
    }
    return 0;
}

