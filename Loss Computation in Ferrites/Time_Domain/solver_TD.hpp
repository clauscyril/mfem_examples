#ifndef SOLVERS_TD_HPP
#define SOLVERS_TD_HPP

#include "../headers/customcurl.hpp"
#include "../headers/Ferrite.hpp"
#include <unordered_set>
#include <mfem.hpp>

using namespace mfem;

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif



void TD_sim(Mesh *mesh, const std::function<real_t(real_t)> &NI_func, real_t t_f, int num_steps, Ferrite ferrite, bool visualization);

void TD_sim_by_flux(Mesh *mesh, const Array<real_t> flux_array, real_t t_f, int num_steps, Ferrite ferrite, bool visualization);

void TD_sim_by_flux(Mesh *mesh, const std::function<real_t(real_t)> &flux_func, real_t t_f, int num_steps, Ferrite ferrite, bool visualization);

void TD_sim_by_fluxH(Mesh *mesh, const Array<real_t> flux_array, real_t t_f, int num_steps, Ferrite ferrite, bool visualization);

void TD_sim_by_fluxH(Mesh *mesh, const std::function<real_t(real_t)> &fluxH_func, real_t t_f, int num_steps, Ferrite ferrite, bool visualization);

void TD_sim_by_v(Mesh *mesh, const std::function<real_t(real_t)> &v_func, real_t t_f, int num_steps, Ferrite ferrite, bool visualization);

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

void readCSV(std::string path, std::vector<std::vector<std::string>> &csvRows);

#endif // SOLVERS_HPP
