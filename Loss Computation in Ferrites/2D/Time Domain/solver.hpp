#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include <mfem.hpp>

using namespace mfem;


// custom curl of an axisymetric problem, with H = H_theta . e_theta and curl(H) = -dH_theta/dz . e_r + (dH_theta/dr + H_tehta/r) e_z
// This class is inspired from GradientGridFunctionCoefficient
class CurlCustomCoefficient : public VectorCoefficient
{
protected:
   const GridFunction *GridFunc;

public:
   CurlCustomCoefficient(const GridFunction *gf);
   void SetGridFunction(const GridFunction *gf);
   const GridFunction * GetGridFunction() const { return GridFunc; }

   void Eval(Vector &V, ElementTransformation &T,
             const IntegrationPoint &ip) override;

   virtual ~CurlCustomCoefficient() { }
};

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


// Fonction qui permet d'imposer les conditions aux limites
real_t bdr_func(const Vector &x);

// Function used for adding the r factor in the integral as we are in cylindrical coordinates
real_t r_coeff_func(const Vector &x);

// Fonction qui permet d'ajouter le terme de compensation en passant en 2D (Passage de rot(rot) à div(grad) + 1/r²)
real_t inv_r_square_func(const Vector &x);


void GetPowerLossByFlux(const char* path, real_t fc, real_t fc_mu, real_t & P_loss_by_vol_mean, std::complex<real_t> &phi);

// Fonction qui calcule la puissance moyenne sur la surface (W/m^3)
void GetPowerLoss(const char* path, real_t fc, real_t fc_mu, real_t &P_loss_eddy, real_t &P_loss_mag, real_t &flux, real_t &Imax);

#endif // SOLVERS_HPP
