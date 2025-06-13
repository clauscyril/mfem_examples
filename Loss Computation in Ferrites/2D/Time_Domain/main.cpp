#include <iostream>
#include <complex>
#include "mfem.hpp"
#include "solver_TD.hpp"

// #include "librom.h"

using namespace mfem;

int main() {
    // const char *path = "../../../mesh/square.msh";
    const char *path = "C:/Users/cyril/Projets/mfem_examples/Loss Computation in Ferrites/mesh/square.msh";
    Mesh *mesh = new Mesh(path, 1, 1); 

    Ferrite N30("N30", 5.98e-2, 4.44e-1, 2.48e-6, 4300);
    Ferrite N87("N87", 4.24e-2, 1.48e-1, 2.68e-6, 2200);
    Ferrite T38("T38", 4.04e-2, 1.07e1, 8.06e-6, 10000);
    
    // real_t I_rms = 0.082348/sqrt(2);
    real_t f = 500e3;
    real_t I_rms = 0.082/sqrt(2);
    int nb_period = 5;
    int nb_iter = 100 * nb_period;
    real_t Ts = nb_period/f/nb_iter;

    real_t t_f = nb_period/f;

    // 
    auto NI_sine_func = [&](real_t t) 
    {
        return I_rms * std::sqrt(2) * std::sin(2 * M_PI * f * t);
    };

    auto NI_square_func = [&](real_t t) 
    {
        int iter = int(t/Ts);
        int rest = iter%((int)(nb_iter/nb_period));

        if (rest < nb_iter/nb_period/2)
            return I_rms * std::sqrt(2)/2;
        else
            return -I_rms * std::sqrt(2)/2.;
    };

    
    auto NI_square_pos_func = [&](real_t t) 
    {
        int iter = int(t/Ts);
        int rest = iter%((int)(nb_iter/nb_period));

        if (rest < nb_iter/nb_period/2)
            return I_rms * std::sqrt(2);
        else
            return 0.;
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
    auto NI_zero = [&](real_t t) 
    {
        return 0;
    };


    TD_sim(mesh, NI_saw_func, t_f, nb_iter, N30, false);
    return 0;
}
