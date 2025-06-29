#include <iostream>
#include <complex>
#include "mfem.hpp"

#include "../headers/customcurl.hpp"
#include "solver.hpp"

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif

using namespace mfem;

// Fonction pour lire un fichier CSV et stocker les données dans un vecteur de doubles
std::vector<std::vector<double>> readCSV(const std::string& filename);


int main(){
    const char *path = "../../../mesh/square2.msh";             // Path to the mesh
    
    // // If generated by mfem :
    // int nr = 10;         // Number of nodes on r axis
    // int nz = 15;         // Number of nodes on y axis

    // // Generation of the mesh
    // Mesh *mesh = new Mesh(nr, nz, Element::TRIANGLE, true, 5.3e-3, 7.59e-3);

    // real_t dx = 4.8e-3;  // r offset

    // for (int i = 0; i < mesh->GetNV(); i++)
    // {
    //     real_t *v = mesh->GetVertex(i);
    //     v[0] += dx;  // addition of the offset for each node
    // }

    Mesh *mesh = new Mesh(path, 1, 1);


    // ********************************
    // // Path to csv files for saving results
    std::string name = "./datafolder/data_tau_fixed.csv";   // Path to csv file for python plot
    std::ofstream data_file(name);                          // ofstream for writing in the file
    data_file << "fc;P_eddy;P_mag;P_tot;flux_r;flux_i;Imax\n";                    // Intialising the file with coluns names

    // Frequency range for the simulation
    real_t fc_0 = 1e3;
    real_t fc = fc_0;
    real_t fc_end = 2e6;
    int N = 49; // N+1 points

    // Change of variable to obtain points linearly spaced on a log scale (optional)
    real_t u = log(fc_0);
    real_t u_end = log(fc_end);
    real_t delta_u = (u_end - u)/ N; 

    real_t fc_mu = 1.8e6;  // cutoff frequency of µeq = µ0 µr / (1 + jw/wc) with wc = 2 pi fc

    for (int i = 0; i < N + 1; i++) {    // Looping the differents frequencies
        fc = exp(u + i*delta_u);                            // Frequency for the simulation (from the change of variable)
        real_t PLoss_eddy, PLoss_mag;
        std::complex<real_t> flux(0,0);                                // parameters to be computed (will be passed as reference)
        real_t imax = 0;  // If not 0, is used for the boundary conditions, else, the flux is set as the condition by rescaling

        GetPowerLoss(mesh, fc, fc_mu, PLoss_eddy, PLoss_mag, flux, imax, false);   // Calling the functions that computes the power losses
        data_file << fc << ";" << PLoss_eddy << ";" << PLoss_mag << ";" << PLoss_eddy + PLoss_mag << ";" << flux.real() << ";" << flux.imag() << ";" << imax << std::endl; // Writing the results in the csv file
    }




    // *****************************************

    // std::vector<std::vector<real_t>> f = readCSV("f.csv");
    // std::vector<std::vector<real_t>> I = readCSV("I.csv");
    // std::string name = "./data";
    // name += std::to_string(0) + ".csv";
    
    // // ofstream to those files
    // std::ofstream data_file(name);
    // data_file << "fc;Ploss;flux\n";

    // for (int i = 0; i<f.size(); i++) {
    //     std::cout << "f = " << f[i][0] << ", I = " << I[i][0] << std::endl;
    //     real_t fc = f[i][0] * 1000;
    //     real_t Imax = I[i][0];
    //     real_t fc_mu = 1e6;

    //     real_t PLoss, flux;

    //     GetPowerLoss(path, fc, fc_mu, PLoss, flux, Imax);
    //     data_file << fc << ";" << PLoss << ";" << flux <<std::endl; 
    // }


    // *****************************************


    // real_t P_eddy, P_mag;
    // std::complex<real_t> flux = 0;
    // real_t imax = 0;
    // GetPowerLoss(mesh, 100e1, 1.8e6, P_eddy, P_mag, flux, imax, true);

return 0;
}



std::vector<std::vector<double>> readCSV(const std::string& filename) {
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier " << filename << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream lineStream(line);
        std::string cell;

        while (std::getline(lineStream, cell, ',')) {
            // Convertir la cellule en double et l'ajouter à la ligne
            row.push_back(std::stod(cell));
        }

        data.push_back(row);
    }

    file.close();
    return data;
}