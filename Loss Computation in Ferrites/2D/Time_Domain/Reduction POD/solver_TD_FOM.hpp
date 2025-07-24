#ifndef SOLVERS_TD_FOM_HPP
#define SOLVERS_TD_FOM_HPP

#include "../../headers/customcurl.hpp"
#include "../../headers/Ferrite.hpp"
#include <mfem.hpp>
#include "linalg/BasisGenerator.h"
#include "linalg/BasisReader.h"
#include <unordered_set>

using namespace mfem;

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif

/**
 * @brief Time-domain electromagnetic simulation with snapshot generation for model order reduction.
 *
 * This function performs a time-domain finite element simulation of the electromagnetic response 
 * of a 2D axisymmetric structure (in cylindrical coordinates) using a parallel mesh `pmesh`.
 * It models the evolution of the magnetic and electric fields in a ferrite material, 
 * applying discrete-time Z-transform formulations and finite element discretization in space.
 *
 * The function computes:
 *  - Magnetic field H and its time evolution.
 *  - Magnetic flux density B.
 *  - Electric field E, derived from the current density J.
 *  - Eddy current losses (power dissipation).
 *
 * It handles:
 *  - Time-dependent Dirichlet boundary conditions driven by the input excitation `NI_func`.
 *  - Optional real-time visualization using GLVis.
 *  - Optional snapshot saving (for reduced-order modeling with libROM).
 *
 * @param pmesh           Parallel finite element mesh of the geometry.
 * @param NI_func         Time-dependent excitation function (e.g., NI(t)).
 * @param t_f             Final simulation time.
 * @param num_steps       Number of time steps.
 * @param ferrite         Material properties (rho, mu, sigma, epsilon).
 * @param visualization   Enable real-time visualization (GLVis) if true.
 * @param save_snapshot   Enable saving field snapshots for model reduction if true.
 * @param num_procs       Total number of MPI processes.
 * @param myid            Rank of the current MPI process.
 */
void TD_sim_offline(ParMesh &pmesh, const std::function<real_t(real_t)> &NI_func, real_t t_f, int num_steps, Ferrite ferrite, bool visualization, bool save_snapshot, int num_procs, int myid);

/**
 * @brief Performs a reduced-order time-domain simulation of electromagnetic behavior in a ferrite component.
 *
 * This function simulates the magnetic and electric fields in a 2D axisymmetric ferrite structure using a reduced-order model (ROM) 
 * based on snapshots obtained offline from the function above. It assembles reduced matrices, projects the problem into a reduced space, 
 * and computes the magnetic field, electric field, current density, and power losses at each time step. 
 * The results can be visualized and are written to a file for post-processing.
 *
 * @param pmesh         Parallel mesh representing the geometry.
 * @param NI_func       Time-dependent source function (NI(t)).
 * @param t_f           Final simulation time.
 * @param num_steps     Number of time steps.
 * @param ferrite       Material parameters of the ferrite (rho, mu, sigma, epsilon).
 * @param visualization Flag to enable/disable visualization.
 * @param num_procs     Number of MPI processes.
 * @param myid          ID of the current MPI process.
 */
void TD_sim_online(ParMesh &pmesh, const std::function<real_t(real_t)> &NI_func, real_t t_f, int num_steps, Ferrite ferrite, bool visualization, int num_procs, int myid);


/*
    Interpolation class for computing the radial component of the curl in cylindrical coordinates
    for axisymmetric problems. In such cases, the vector field is assumed to have only a theta
    component, and all derivatives with respect to theta are zero.

    This class interpolates the radial component of the curl using only the theta component of the field.
    The interpolation is performed in an H1 nodal space, which may introduce errors since the true
    curl is generally discontinuous.

    For each degree of freedom, the interpolated value is taken from the first element that contributes
    to it. As a result, the interpolation may exhibit artifacts or noise near element boundaries. 
*/
class CustomCurlInterpolatorX : public DiscreteInterpolator
{
public:
   CustomCurlInterpolatorX() {};
   virtual ~CustomCurlInterpolatorX() {};


/* 
 * Computes the local interpolation matrix on the degrees of freedom (DOFs) 
 * associated with a given element.
 *
 * This function is invoked by ParDiscreteLinearOperator::Assemble(), 
 * which is inherited from BilinearForm. During the global assembly process, 
 * Assemble() iterates over all mesh elements, computes the corresponding local 
 * matrices via this function, and assembles them into the global matrix.
 *
 * Parameters:
 *  - h1_fe:    The finite element representing the H1 space on the current element.
 *  - out_fe:   The output finite element space used for interpolation.
 *  - Trans:    The element transformation providing geometric mappings.
 *  - elmat:    The resulting local interpolation matrix to be filled.
 */
void AssembleElementMatrix2(const FiniteElement &h1_fe,
                            const FiniteElement &out_fe,
                            ElementTransformation &Trans,
                            DenseMatrix &elmat) override;

private:
    std::unordered_set<real_t> set;  // Set used to identify the different degrees of freedom 
};

// The Class does exactly the same as the above one, only interpolating the z componant of the curl
class CustomCurlInterpolatorY : public DiscreteInterpolator
{
public:
   CustomCurlInterpolatorY() {};
   virtual ~CustomCurlInterpolatorY() {};

void AssembleElementMatrix2(const FiniteElement &h1_fe,
                            const FiniteElement &out_fe,
                            ElementTransformation &Trans,
                            DenseMatrix &elmat) override;

private:
    std::unordered_set<real_t> set;
};

// Coefficient for evaluating the losses in Time Domain. 
class PowerLossCoefficient_TD : public mfem::Coefficient
{
private:
    const FiniteElementSpace *fespace;
    CurlCustomCoefficient J;
    mutable mfem::Vector J_vect;
    VectorGridFunctionCoefficient E;
    mutable mfem::Vector E_vect;

public:
    PowerLossCoefficient_TD(const FiniteElementSpace *fespace_, CurlCustomCoefficient &J_, VectorGridFunctionCoefficient &E_);
    virtual real_t Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override;
};

// Compute the Matrix product C^t A B, with A a mfem::SparseMatrix
// Inspired from ComputeCtAB of librom : https://librom.readthedocs.io/en/latest/ for mfem::HypreParMatrix
void ComputeCtAB_Sparse(const mfem::SparseMatrix& A,
                        const CAROM::Matrix& B,  // Distributed matrix
                        const CAROM::Matrix& C,  // Distributed matrix
                        CAROM::Matrix& CtAB);     // Non-distributed (local) matrix

#endif // SOLVERS_HPP
