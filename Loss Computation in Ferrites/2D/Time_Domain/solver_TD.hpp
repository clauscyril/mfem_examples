#ifndef SOLVERS_TD_HPP
#define SOLVERS_TD_HPP

#include "../headers/customcurl.hpp"
#include "../headers/Ferrite.hpp"
#include <mfem.hpp>

using namespace mfem;

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif



void TD_sim(Mesh *mesh, const std::function<real_t(real_t)> &NI_func, real_t t_f, int num_steps, Ferrite ferrite, bool visualization);


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
