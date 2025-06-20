#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include <mfem.hpp>
#include "../headers/customcurl.hpp"
#include "../headers/Ferrite.hpp"

using namespace mfem;


// Class personnalisée permettant de calculer la densité de puissance Re(E.J*) = Re(rho) (||Jr||² + ||Ji||²) 
class PowerLossCoefficient : public mfem::Coefficient
{
private:
    const FiniteElementSpace *fespace;
    CurlCustomCoefficient J_r;
    mutable mfem::Vector J_r_vect;
    CurlCustomCoefficient J_i;
    mutable mfem::Vector J_i_vect;
    std::complex<real_t> rho_eq;

public:
    PowerLossCoefficient(const FiniteElementSpace *fespace_, std::complex<real_t> rho_eq_, CurlCustomCoefficient &J_r_, CurlCustomCoefficient &J_i_);
    virtual real_t Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override;
};

// Class personnalisée permettant de calculer la densité de puissance Re(E.J*) = Re(rho) (||Jr||² + ||Ji||²) 
class PowerLossMagCoefficient : public mfem::Coefficient
{
private:
    const FiniteElementSpace *fespace;
    GridFunctionCoefficient H_r;
    GridFunctionCoefficient H_i;
    std::complex<real_t> mu_eq;
    real_t omega;

public:
    PowerLossMagCoefficient(const FiniteElementSpace *fespace_, std::complex<real_t> mu_eq_, real_t omega_, GridFunctionCoefficient &H_r_, GridFunctionCoefficient &H_i_);
    virtual real_t Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override;
};



// Fonction qui calcule la puissance moyenne sur la surface (W/m^3)
void GetPowerLoss(Mesh* mesh, real_t fc, real_t fc_mu, real_t &P_loss_eddy, real_t &P_loss_mag, std::complex<real_t> &flux, real_t &Imax, const bool visualization);

#endif // SOLVERS_HPP
