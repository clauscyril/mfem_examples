#include <iostream>
#include <filesystem>
#include <complex>
#include "mfem.hpp"
#include <vector>

#include "../headers/customcurl.hpp"
#include "solver.hpp"

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif

using namespace mfem;

void readCSV(std::string path, std::vector<std::vector<std::string>> &csvRows) {
  	std::ifstream input{path};

	if (!input.is_open()) {
		std::cerr << "Couldn't read file: " << path << "\n";
		return; 
	}

	for (std::string line; std::getline(input, line);) {
		std::istringstream ss(std::move(line));
		std::vector<std::string> row;
		if (!csvRows.empty()) {
			// We expect each row to be as big as the first row
			row.reserve(csvRows.front().size());
		}

		for (std::string value; std::getline(ss, value, ',');) {
			row.push_back(std::move(value));
		}
		csvRows.push_back(std::move(row));
	}
	return;
}

int main(int argc, char *argv[]){        
    const char *path = "../../mesh/N30_P.msh";
    Mesh *mesh = new Mesh(path, 1,1);
    mesh->UniformRefinement();

    std::string material = "N30";
    Ferrite N30("N30", 5.98e-2, 4.44e-1, 2.48e-6, 4832, 1200e3);
    Ferrite N87("N87", 4.24e-2, 1.48e-1, 2.68e-6, 2200, 4e6);
    Ferrite T38("T38", 4.04e-2, 1.07e1, 8.06e-6, 10000, 380e3);
    
    Ferrite ferrite = N30;

    
    std::string path_data = "../data/Princeton/N30/N30-Sinusoidal_Phi_15.csv";
    OptionsParser args(argc, argv);
    args.AddOption(&material, "-m", "--material", "Material to use.");
    args.AddOption(&path_data, "-p", "--path_data", "path to measurements data");
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

    std::filesystem::path fs_path = path_data;
    if (!std::filesystem::exists(fs_path)) {
        "File not found, using data/Princeton/N30/N30-Sinusoidal_Phi_15.csv instead";
        path_data = "../data/Princeton/N30/N30-Sinusoidal_Phi_15.csv";
    }

    // Rows from csv file
    std::vector<std::vector<std::string>> Rows;
    readCSV(path_data, Rows);

    std::string phi = std::to_string(static_cast<int>(1000*std::stod(Rows[1][1])));  // Get b_peak value in mT in a string
    std::string name = "../data/MFEM/Princeton_mfem"+ material + "_"+ phi + ".csv"; 

    std::ofstream data_file(name);                         
    data_file << "fc;P_eddy;P_mag;P_tot;flux_r;flux_i\n";

    for (int i = 1; i<Rows.size(); i++){
        real_t f =  std::stod(Rows[i][0]);
        real_t b_peak = std::stod(Rows[i][1]);   

        real_t PLoss_eddy, PLoss_mag;
        std::complex<real_t> flux(0,0);  
        real_t NI = 0;

        GetPowerLoss(mesh, f, b_peak, NI, ferrite, PLoss_eddy, PLoss_mag, flux, false); 
        data_file << f << ";" << PLoss_eddy << ";" << PLoss_mag << ";" << PLoss_eddy + PLoss_mag << ";" << flux.real() << ";" << flux.imag()  << std::endl;
    }   


return 0;
}