#include "mfem.hpp"
#include "librom.h"
#include <fstream>
#include <iostream>
#include "linalg/BasisGenerator.h"
#include "linalg/BasisReader.h"
#include "mfem/Utilities.hpp"

#include "solver_TD_FOM.hpp"

using namespace std;
using namespace mfem;

// rhs function. See below for implementation.

int dim;

int main(int argc, char *argv[])
{

    Ferrite N30("N30", 5.98e-2, 4.44e-1, 2.48e-6, 4300);
    Ferrite N87("N87", 4.24e-2, 1.48e-1, 2.68e-6, 2200);
    Ferrite T38("T38", 4.04e-2, 1.07e1, 8.06e-6, 10000);

    const char *mesh_file = "../../../../mesh/square.msh";

    Mesh mesh(mesh_file,1,1);
    mesh.UniformRefinement();
    // mesh.UniformRefinement();
    // mesh.UniformRefinement();

    // Fonction source (Courant dans la bobine sous forme de condition aux limites)
    real_t f = 800e3;
    real_t I_rms = 0.082/sqrt(2);
    int nb_period = 3;
    int nb_iter = 100 * nb_period;
    real_t Ts = nb_period/f/nb_iter;

    real_t t_f = nb_period/f;
    auto NI_sine_func = [&](real_t t) 
    {
        return I_rms * std::sqrt(2) * std::sin(2 * M_PI * f * t);
    };

    auto NI_saw_pos_func = [&](real_t t) 
    {
        int iter = int(t/Ts);
        int rest = iter%((int)(nb_iter/nb_period));

        return I_rms * std::sqrt(2) * rest * Ts;
    };  

    auto NI_saw_func = [&](real_t t) 
    {
        int iter = int(t/Ts);
        int rest = iter%((int)(nb_iter/nb_period));

        if (rest < nb_iter/nb_period/2)
            return I_rms * std::sqrt(2) * rest * Ts;
        else
            return I_rms * std::sqrt(2) * (nb_iter/nb_period - rest) * Ts;
    };

    // TD_sim_offline(mesh, NI_sine_func, t_f, nb_iter, N30, false, true);
    // std::cout << "test" << std::endl;
    TD_sim_online(mesh, NI_sine_func, t_f, nb_iter, N30, false);

    return 0;
}

