#include <iostream>
#include <complex>
#include "mfem.hpp"
#include "solver_TD.hpp"

using namespace mfem;

int main() {
    const char *path = "../../../mesh/square.msh";
    Mesh *mesh = new Mesh(path, 1, 1); 

    Ferrite N30("N30", 5.98e-2, 4.44e-1, 2.48e-6, 4300);
    Ferrite N87("N87", 4.24e-2, 1.48e-1, 2.68e-6, 2200);
    Ferrite T38("T38", 4.04e-2, 1.07e1, 8.06e-6, 10000);
    
    real_t I_rms = 0.082348/sqrt(2);
    real_t f = 100e3;
    real_t Ts = 1/f / 100;

    TD_sim(mesh, I_rms, f, Ts, 1000, N30, false);
    return 0;
}
