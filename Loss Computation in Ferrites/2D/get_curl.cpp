#include "get_curl.hpp"
#include <iostream>
#include <cmath>


using namespace mfem;

CurlCustomCoefficient::CurlCustomCoefficient (
   const GridFunction *gf)
   : VectorCoefficient((gf) ?
                       gf -> FESpace() -> GetMesh() -> SpaceDimension() : 0)
{
   GridFunc = gf;
}

void CurlCustomCoefficient::SetGridFunction(const GridFunction *gf)
{
   GridFunc = gf; vdim = (gf) ?
                         gf -> FESpace() -> GetMesh() -> SpaceDimension() : 0;
}

void CurlCustomCoefficient::Eval(Vector &V, ElementTransformation &T,
                                           const IntegrationPoint &ip)
{
    Vector grad_H;
    GridFunc->GetGradient(T, grad_H);
    V.SetSize(2);
    V[0] = grad_H[1];
    V[1] = -grad_H[0];
}




PowerLossCoefficient::PowerLossCoefficient(const FiniteElementSpace *fespace_, std::complex<real_t> rho_eq_, CurlCustomCoefficient &J_r_, CurlCustomCoefficient &J_i_)
    : fespace(fespace_),
      J_r(J_r_), J_r_vect(fespace_->GetMesh()->SpaceDimension()),
      J_i(J_i_), J_i_vect(fespace_->GetMesh()->SpaceDimension()),
      rho_eq(rho_eq_) {}

real_t PowerLossCoefficient::Eval(ElementTransformation &T,
                                    const IntegrationPoint &ip)
{
    J_r.Eval(J_r_vect, T, ip);
    J_i.Eval(J_r_vect, T, ip);
    return rho_eq.real() * (pow(J_r_vect.Norml2(),2)+ pow(J_i_vect.Norml2(),2));
}