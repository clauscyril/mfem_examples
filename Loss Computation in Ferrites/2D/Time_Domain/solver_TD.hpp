#ifndef SOLVERS_TD_HPP
#define SOLVERS_TD_HPP

#include "../customcurl.hpp"
#include <mfem.hpp>

using namespace mfem;

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif



struct Ferrite
{
    std::string name;
    real_t rho;
    real_t sigma;
    real_t eps;
    real_t mu; 
    Ferrite(const char* name_, real_t rho_, real_t sigma_, real_t eps_, real_t mu_r_); 
};

void TD_sim(Mesh *mesh, real_t I_rms, real_t f, real_t Ts, int num_steps, Ferrite ferrite, bool visualization);


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

#endif // SOLVERS_HPP
