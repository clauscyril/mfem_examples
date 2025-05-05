#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include <mfem.hpp>

using namespace mfem;

/// Fonction source temporaire
mfem::real_t q_func(const mfem::Vector &x, mfem::real_t &t);

/// Solver implicite DG
class DG_Solver : public Solver
{
private:
   SparseMatrix &M, &K, A;
   GMRESSolver linear_solver;
   BlockILU prec;
   real_t dt;
public:
   DG_Solver(SparseMatrix &M_, SparseMatrix &K_, const FiniteElementSpace &fes);
   void SetTimeStep(real_t dt_);
   void SetOperator(const Operator &op) override;
   void Mult(const Vector &x, Vector &y) const override;
};

/// Opérateur d’évolution pour le FEM
class FE_Evolution : public TimeDependentOperator
{
private:
   BilinearForm &M, &K;
   const Vector &b;
   Solver *M_prec;
   CGSolver M_solver;
   DG_Solver *dg_solver;
   mutable Vector z;
public:
   FE_Evolution(BilinearForm &M_, BilinearForm &K_, const Vector &b_);
   void Mult(const Vector &x, Vector &y) const override;
   void ImplicitSolve(const real_t dt, const Vector &x, Vector &k) override;
   ~FE_Evolution() override;
};

#endif // SOLVERS_HPP
