#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include <mfem.hpp>

using namespace mfem;

class CurlCustomCoefficient : public VectorCoefficient
{
protected:
   const GridFunction *GridFunc;

public:

   /** @brief Construct the coefficient with a scalar grid function @a gf. The
       grid function is not owned by the coefficient. */
   CurlCustomCoefficient(const GridFunction *gf);

   ///Set the scalar grid function.
   void SetGridFunction(const GridFunction *gf);

   ///Get the scalar grid function.
   const GridFunction * GetGridFunction() const { return GridFunc; }

   /// Evaluate the gradient vector coefficient at @a ip.
   void Eval(Vector &V, ElementTransformation &T,
             const IntegrationPoint &ip) override;

   /** @brief Evaluate the gradient vector coefficient at all of the locations
       in the integration rule and write the vectors into columns of matrix @a
       M. */
//    void Eval(DenseMatrix &M, ElementTransformation &T,
//              const IntegrationRule &ir) override;

   virtual ~CurlCustomCoefficient() { }
};


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
    virtual real_t Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
};

#endif // SOLVERS_HPP



real_t bdr_func(const Vector &x);
real_t inv_r_square_func(const Vector &x);


real_t GetPowerLoss(const char *path, real_t fc, real_t fc_mu);

