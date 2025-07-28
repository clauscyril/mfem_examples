#include "mfem.hpp"
#include "librom.h"
#include <fstream>
#include <iostream>
#include "linalg/BasisGenerator.h"
#include "linalg/BasisReader.h"
#include "mfem/Utilities.hpp"

#include "solver_TD_FOM.hpp"

using namespace mfem;

int main(int argc, char *argv[])
{

    Ferrite N30("N30", 5.98e-2, 4.44e-1, 2.48e-6, 4300, 1.8e6);
    Ferrite N87("N87", 4.24e-2, 1.48e-1, 2.68e-6, 2200, 1.8e6);
    Ferrite T38("T38", 4.04e-2, 1.07e1, 8.06e-6, 10000, 1.8e6);

    const char *mesh_file = "../../../../mesh/square.msh";

    Mesh mesh(mesh_file,1,1);
    mesh.UniformRefinement();
    mesh.UniformRefinement();

    // Time sumlation settings
    real_t f = 100e3;                   // Frequency of the signal 
    real_t I_rms = 0.082/sqrt(2);       // RMS Value
    int nb_period = 5;                  // Number of period simulated
    int nb_iter = 100 * nb_period;      // 100 time iterations for each period
    real_t Ts = nb_period/f/nb_iter;    
    real_t t_f = nb_period/f;

    // -------- Different source functions for the current------------

    // Simple sine wave
    auto NI_sine_func = [&](real_t t) 
    {
        return I_rms * std::sqrt(2) * std::sin(2 * M_PI * f * t);
    };

    // Saw wave form
    auto NI_saw_func = [&](real_t t) 
    {
        int iter = int(t/Ts);
        int rest = iter%((int)(nb_iter/nb_period));
        return I_rms * std::sqrt(2) * rest * Ts;
    };  

    // Triangle waveform
    auto NI_triangle_func = [&](real_t t) 
    {
        int iter = int((t + 0.25/f)/Ts);
        int rest = iter%((int)(nb_iter/nb_period));
        real_t Ipeak = I_rms * sqrt(3);
        if (rest < nb_iter/nb_period/2)
            return 2 * Ipeak * rest/(nb_iter/nb_period/2) - Ipeak;
        else
            return - 2 * Ipeak * rest/(nb_iter/nb_period/2) + 3* Ipeak;
    };
    // ----------------------------------------------------------------

    // ---------- Initializing MPI -----------
    Mpi::Init();
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    Hypre::Init();

    ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();
    // ---------------------------------------
    

    // Offline simulation : Full FEM simulation, used for generating the snapshots for the reduced model 
    TD_sim_offline(pmesh, NI_saw_func, t_f, nb_iter, N30, false, false, num_procs, myid);
    
    // Online simulation : Uses the snaphots generated to generate the reduced space.
    TD_sim_online(pmesh, NI_saw_func, t_f, nb_iter, N30, false, num_procs, myid);
    
    Mpi::Finalize();
    return 0;
}

