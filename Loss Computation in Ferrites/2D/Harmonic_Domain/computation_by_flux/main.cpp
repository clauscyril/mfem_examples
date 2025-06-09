#include <iostream>
#include <complex>
#include "mfem.hpp"

#include "solver.hpp"

#define M_PI 3.14159265358979323846 

using namespace mfem;

// Fonction pour lire un fichier CSV et stocker les données dans un vecteur de doubles
std::vector<std::vector<double>> readCSV(const std::string& filename);


int main(){
    const char *path = "../../mesh/square.msh";             // Path to the mesh


    // *****************************************

    std::vector<std::vector<real_t>> f = readCSV("f.csv");
    std::vector<std::vector<real_t>> I = readCSV("I.csv");
    std::string name = "./data";
    name += std::to_string(0) + ".csv";
    
    // ofstream to those files
    std::ofstream data_file(name);
    data_file << "fc;Ploss;flux\n";

    for (int i = 0; i<f.size(); i++) {
        std::cout << "f = " << f[i][0] << ", I = " << I[i][0] << std::endl;
        real_t fc = f[i][0] * 1000;
        real_t Imax = I[i][0];
        real_t fc_mu = 1e6;

        real_t PLoss, flux;

        GetPowerLoss(path, fc, fc_mu, PLoss, flux, Imax);
        data_file << fc << ";" << PLoss << ";" << flux <<std::endl; 
    }


    // *****************************************


    real_t PLoss;
    std::complex<real_t> phi(1,1);
    
    // GetPowerLossByFlux(path, 0.5e6, 3e6, PLoss, phi);
    // std::cout << PLoss << std::endl;

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