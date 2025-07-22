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



void TD_sim_offline(Mesh &mesh, const std::function<real_t(real_t)> &NI_func, real_t t_f, int num_steps, Ferrite ferrite, bool visualization, bool save_snapshot);

void TD_sim_online(Mesh &mesh, const std::function<real_t(real_t)> &NI_func, real_t t_f, int num_steps, Ferrite ferrite, bool visualization);


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
